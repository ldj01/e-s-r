#!/usr/bin/env python3

import sys
import os
import logging
from datetime import datetime, timedelta

import cx_Oracle
from espa import System
from espa import Date
from espa import FtpSession
from espa import AuxConfig

# Global variables
ERROR = 1
SUCCESS = 0
logger = logging.getLogger(__name__)

def getNcepData(cfg, year, overwrite):
    """ Downloads the air temp, surface pressure, and precipitable water netCDF
        files for the desired year, then process the netCDF file into individual
        daily HDF files containing the air temp, surface pressure, and
        precipitable water.

        Parameters:
            cfg: config information for filename formats
            year: year of NCEP data to be downloaded and processed (integer)
            overwrite: boolean to determine whether existing archive files are
                       overwritten or skipped

        Returns:
            SUCCESS/ERROR on processing success or failure
    """
    # gather directory and filename format info
    ancdir = cfg.get('ncep_base_archive')
    dloaddir = cfg.get('ncep_temp_directory')
    outputDir = cfg.get('ncep_archive_format').format(ancdir, year)
    archiveFileFormat = cfg.get('ncep_filename_format')

    # set up the URLs of the NCEP netCDF files to be downloaded for the
    # specified year [pressure file, water file, air file]
    netCDF_list = [cfg.get('pressure_url_format').format(year),
                   cfg.get('water_url_format').format(year),
                   cfg.get('air_url_format').format(year)]

    # connect to the database
    try:
        con = cx_Oracle.connect(os.getenv('IAS_DB_COM'))
        dbh = con.cursor()
        logger.info('Connected to database successfully')
    except cx_Oracle.DatabaseError as e:
        logger.error('Could not connect to database: %s', e)
        return ERROR

    # make sure the path exists and is cleaned up
    System.create_directory(dloaddir)
    System.empty_directory(dloaddir)

    # download the air temp, surface pressure, and precipitable water files
    # for the specified year to the temporary download directory
    session = FtpSession()
    for netCDF in netCDF_list:
        status = session.ftp_transfer_file(netCDF, dloaddir)
        if status != SUCCESS:
            logger.error('could not download file: %s', netCDF)
            return ERROR

    # use the downloaded netCDF files to create the daily HDF files needed
    # for LEDAPS processing
    # make sure the output directory exists or create it recursively
    if not os.path.exists(outputDir):
        logger.warn('%s does not exist... creating', outputDir)
        System.create_directory(outputDir)

    # if the specified year is the current year, only process up through
    # today otherwise process through all the days in the year
    day_of_year = Date.get_days_elapsed_in_year(year)

    # generate a list of the current archive
    archive = System.directory_file_list(outputDir)

    # loop through each day in the year and process the NCEP REANALYSIS HDF
    # file for each day
    for doy in range(1, day_of_year + 1):
        fname = archiveFileFormat.format(year, doy)
        # skip files that do not need to be overwritten
        if not overwrite and os.path.join(outputDir, fname) in archive:
            continue

        for netCDF in netCDF_list:
            netCDF_name = os.path.join(dloaddir, os.path.basename(netCDF))
            status = executeNcep(netCDF_name, os.path.join(outputDir, fname),
                                 year, doy)
            if status == ERROR:
                logger.error('could not process file: %s', netCDF_name)
                return ERROR
        # build the db entry for this file and insert/update the record
        values = {}
        values['coverage'] = (datetime(year, 1, 1) + timedelta(days=doy-1))
        values['fname'] = fname

        try:
            existing = overwrite and os.path.join(outputDir, fname) in archive
            if existing:
                # if we are overwriting and the file already exists, then a db
                # record already exists and must be updated
                logger.info('Updating existing database archive entry')
                dbh.execute('UPDATE REANALYSIS SET FILE_NAME = :fname, '
                            'DATE_ENTERED = SYSDATE WHERE EFFECTIVE_DATE = '
                            ':coverage', values)
                if dbh.rowcount != 1:
                    logger.error('No rows updated')
                    raise cx_Oracle.DatabaseError
            else:
                logger.info('Creating new database archive entry')
                dbh.execute('INSERT INTO REANALYSIS (EFFECTIVE_DATE, '
                            'FILE_NAME, DATE_ENTERED) VALUES (:coverage, '
                            ':fname, SYSDATE)', values)
            con.commit()
            logger.info('Processing successful for DOY %s, year %s', doy, year)
        except cx_Oracle.DatabaseError as e:
            logger.error('Database entry for DOY %s, year %s has failed: %s'
                         'Exiting processing.', doy, year, e)
            # if the file was not already in the archive, remove it
            if not existing:
                logger.error('Removing the archived file to match db records')
                os.remove(os.path.join(outputDir, fname))
            else:
                logger.error('This file previously existed in the local '
                             'archive. Double check that the archive and '
                             'database record for this file match.')

            return ERROR
    # end for doy loop

    # close db connection
    dbh.close()
    con.close()

    # cleanup the downloaded annual netCDF files
    System.empty_directory(dloaddir)

    return SUCCESS


def executeNcep(inputfile, outputfile, year, doy):
    """ Runs the 'ncep' executable to produce the HDF files for the specified
        year which are written to outputdir.  If the specified year is the
        current year, then the days processed will only be up through today.
        If the outputdir directory does not exist, then it is made before
        downloading.

        Parameters:
            inputfile: full path and filename of the NCEP REANALYSIS file for
                       the specified year
            outputfile: the full path and filename format for the final
                        destination of the processed REANALYSIS data
            year: year of NCEP data to be processed
            doy: day of year of NCEP data to be processed

        Returns:
            SUCCESS/ERROR on success or failure of processing the NCEP files
    """
    # prepare the ncep_repackage command
    cmdstr = 'ncep_repackage {} {} {}'.format(inputfile, outputfile, doy)
    status = System.execute_cmd(cmdstr)
    if status == 157:  # return value of -99 (2s complement of 157)
        logger.error('Input file for year %s, DOY %s is not'
                     ' readable.  Stop processing since this same file is'
                     ' used for all days in the current year.', year, doy)
        return ERROR
    elif status != 0:
        logger.warn('error running ncep for year %s, DOY %s.  '
                    'Processing will continue ...', year, doy)
        # if processing fails, then remove the output file.
        if os.path.isfile(outputfile):
            os.remove(outputfile)

    return SUCCESS

def main():
    """ Main routine which grabs the command-line arguments, determines
        which dates and times of data need to be processed, then processes the
        NCEP data that falls between the specified dates

        Returns:
            SUCCESS/ERROR on success or failure of MERRA-2 data retrieval and
            processing
    """
    # get ncep config information
    c = AuxConfig()
    cfg = c.get_config('ncep')
    start_date = datetime.strptime(cfg.get('ncep_start_date'), '%Y-%m-%d')

    # gather command line arguments
    description = ('Downloads NCEP data and extracts the required parameters '
                   'for Surface Reflectance processing.  The parameters '
                   'are then archived for later use.')
    try:
        args = System.get_command_line_arguments(description, start_date, True)
    except Exception as e:
        print('Argument error: ', e)
        return ERROR


    # setup the default logger format and level. log to STDOUT.
    logging_level = os.getenv('IAS_LOG_LEVEL')
    if logging_level is None:
        logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging_level)

    # Alert user if the overwrite flag is set
    if args.overwrite:
        logger.info('overwrite flag is set: Any existing archive data will '
                    'be overwritten by most current online archive data.')

    # Iterate through the given date range and process data
    logger.info('Processing NCEP data for %s - %s', args.start_year,
                args.end_year)
    for year in range(args.start_year, args.end_year+1):
        logger.info('Processing year: %s', year)
        status = getNcepData(cfg, year, args.overwrite)
        if status == ERROR:
            logger.warn('Problems occurred while processing NCEP'
                        ' data for year %s.  Processing will continue.', year)

    logger.info('NCEP processing complete.')
    return SUCCESS

if __name__ == "__main__":
    sys.exit(main())
