#!/usr/bin/env python3

import sys
import os
import fnmatch
import re
import logging
from datetime import datetime, timedelta
import cx_Oracle

from espa import System
from espa import Date
from espa import FtpSession
from espa import AuxConfig

# global variables
logger = logging.getLogger(__name__)
ERROR = 1
SUCCESS = 0

def resolve(cfg, year):
    """ Resolve which instrument will be used for the specified year. Identify
        the appropriate URL for downloading data for that instrument, and put
        that url on the list.

        Parameters:
            year: year of desired ozone data

        Returns:
            None: error resolving the instrument and associated URL for
            the specified year
            urlList: List of URLs to pull the ozone data from for the specified
                     year.  The primary data
    """
    urlList = []     # create empty data source list

    # use NIMBUS data for 1978-1990
    if year in range(1978, 1991):
        urlList.append(cfg.get('nimbus_url_format'))

    # use METEOR3 data for 1991-1993, with NIMBUS as the backup
    elif year in range(1991, 1994):
        urlList.append(cfg.get('meteor3_url_format'))
        urlList.append(cfg.get('nimbus_url_format'))

    # use METEOR3 data for 1994
    elif year == 1994:
        urlList.append(cfg.get('meteor3_url_format'))

    # use EARTHPROBE data for 1996-2003
    elif year in range(1996, 2004):
        urlList.append(cfg.get('earthprobe_url_format'))

    # use OMI data for 2004-2005, with EARTHPROBE as the backup
    elif year in range(2004, 2006):
        urlList.append(cfg.get('omi_url_format'))
        urlList.append(cfg.get('earthprobe_url_format'))

    # use OMI for any years beyond 2006
    elif year >= 2006:
        urlList.append(cfg.get('omi_url_format'))

    # year requested does not have TOMS/EP ozone data
    if not urlList or all(url is None for url in urlList):
        logger.warn('Could not resolve a datasource for year: %s', year)
        return None

    return urlList


def resolveFile(fileList):
    """ Looks at the list of available files and grab the highest priority
        file to be used for processing in the order of OMI, EARTHPROBE,
        METEOR3, and NIMBUS7 and return the source of that file.

        Parameters:
            fileList: List of available ozone files for the current day and year

        Returns (filename, ozoneSource):
            filename: Priority file to be processed
            ozoneSource: Source of the priority file to be processed
            none: None of the files matched our known instruments
    """
    omiregex = re.compile(r'.*_omi_\d*.txt')
    earthproberegex = re.compile(r'.*_epc_\d*.txt')
    meteor3regex = re.compile(r'.*_m3t_\d*.txt')
    nimbusregex = re.compile(r'.*_n7t_\d*.txt')

    # loop through the files, looping for OMI, EARTHPROBE, METEOR3, and NIMBUS7
    # files in that order.  return the first one found as the file to be
    # processed.
    for myfile in fileList:
        if omiregex.search(myfile):
            return myfile, 'OMI'

    for myfile in fileList:
        if earthproberegex.search(myfile):
            return myfile, 'EARTHPROBE'

    for myfile in fileList:
        if meteor3regex.search(myfile):
            return myfile, 'METEOR3'

    for myfile in fileList:
        if nimbusregex.search(myfile):
            return myfile, 'NIMBUS7'

    # if none of the files match our known instruments then return None
    return None, None


def downloadToms(cfg, year, datestr, destination):
    """ Retrieves the files for the specified year from the EP/TOMS ftp site
        and download to the desired destination.  If the destination directory
        does not exist, then it is made before downloading. Existing files in
        the download directory are removed/cleaned.

        Parameters:
        year: Year of data to download
        datestr: MMDD format for use in source file download
        destination: Name of the directory on the local system to download the
                     EP/TOMS files
        fname: Filename of the downloaded file on the local system

        Returns:
            SUCCESS/ERROR on successful or failed processing
    """
    # make sure the download directory exists and is empty
    System.create_directory(destination)
    System.empty_directory(destination)

    # obtain the list of URL(s) for our particular date
    urlList = resolve(cfg, year)
    if urlList is None:
        logger.warn('EP/TOMS URL could not be resolved for year %s.'
                    '  Processing will continue ...', year)
        return ERROR

    # create the FtpSession object and download the files to the temp location
    session = FtpSession()
    for url in urlList:
        logger.info('Downloading from %s to %s',
                    url.format(year, year, datestr), destination)
        session.ftp_transfer_file(url.format(year, year, datestr), destination)

    return SUCCESS


def getTomsData(cfg, year, overwrite):
    """ Downloads the daily ozone data files for the desired year, then
        processes the text files into individual daily HDF files containing
        the ozone. If the overwrite flag is not set, any days already in
        the archive will be skipped for downloading.

        Parameters:
        cfg: Config information for filename formats
        year: Year of EP/TOMS data to be downloaded and processed
        overwrite: Boolean to determine whether existing archive files are
                   overwritten or skipped

        Returns:
            SUCCESS/ERROR on successful or failed processing
    """
    # gather directory and filename format info
    ancdir = cfg.get('toms_base_archive')
    dloaddir = cfg.get('toms_temp_directory')
    outputDir = cfg.get('toms_archive_format').format(ancdir, year)
    archiveFileFormat = cfg.get('toms_filename_format')

    # connect to the database
    try:
        con = cx_Oracle.connect(os.getenv('IAS_DB_COM'))
        dbh = con.cursor()
        logger.info('Connected to database successfully')
    except cx_Oracle.DatabaseError as e:
        logger.error('Could not connect to database: %s', e)
        return ERROR

    # if the specified year is the current year, only process up through
    # today otherwise process through all the days in the year
    day_of_year = Date.get_days_elapsed_in_year(year)

    # determine the directory for the output ancillary data files to be
    # processed.  create the directory if it doesn't exist.
    if not os.path.exists(outputDir):
        logger.info('%s does not exist... creating', outputDir)
        System.create_directory(outputDir)

    # generate a list of files in the current archive
    archive = System.directory_file_list(outputDir)

    # loop through each day in the year and process the EP/TOMS data
    for doy in range(1, day_of_year + 1):
        # generate the file and date string associated with the current doy
        fname = archiveFileFormat.format(year, doy)
        datestr = Date.doy_format(year, doy, '%m%d')

        # skip files that do not need to be overwritten
        if not overwrite and os.path.join(outputDir, fname) in archive:
            continue

        # download the daily ozone files for the specified day and year
        status = downloadToms(cfg, year, datestr, dloaddir)
        if status == ERROR:
            return ERROR

        # find all the files for the current day
        fileList = []    # create empty list to store files matching date
        for myfile in os.listdir(dloaddir):
            if fnmatch.fnmatch(myfile, '*' + datestr + ".txt"):
                fileList.append(myfile)

        # make sure files were found or print a warning
        if not fileList: # Empty lists evaluate to False in Python
            logger.warn('no EP/TOMS new data available for doy %s year'
                        ' %s. processing will continue ...', doy, year)
            continue

        # determine the highest priority file found based on source, and the
        # source of that file
        tomsfile, ozoneSource = resolveFile(fileList)
        if tomsfile is None or ozoneSource is None:
            logger.warn('error resolving the list of EP/TOMS'
                        ' files to process. processing will continue ...')
            continue

        # generate the path for the input and output files to be processed
        fullOutputPath = os.path.join(outputDir, fname)
        fullInputPath = os.path.join(dloaddir, tomsfile)

        # prepare and execute the convert_ozone command
        cmd = 'convert_ozone {} {} {}'.format(fullInputPath,
                                              fullOutputPath, ozoneSource)
        status = System.execute_cmd(cmd)
        if status != 0:
            logger.warn("Error executing command '%s': %s", str(cmd), e)

        # build the db entry for this file and insert/update the record
        values = {}
        values['coverage'] = (datetime(year, 1, 1) + timedelta(days=doy-1))
        values['fname'] = fname
        values['instrument'] = ozoneSource

        try:
            existing = overwrite and os.path.join(outputDir, fname) in archive
            if existing:
                # if overwriting and the file is already in the archive,
                # then a db record already exists and must be updated
                logger.info('Updating existing database archive entry')
                dbh.execute('UPDATE EP_TOMS_OZONE SET FILE_NAME = :fname, '
                            'DATE_ENTERED = SYSDATE, INSTRUMENT_SOURCE = '
                            ':instrument WHERE EFFECTIVE_DATE = :coverage',
                            values)
                if dbh.rowcount != 1:
                    logger.error('No rows updated')
                    raise cx_Oracle.DatabaseError
            else:
                logger.info('Creating new database archive entry')
                dbh.execute('INSERT INTO EP_TOMS_OZONE (EFFECTIVE_DATE, '
                            'FILE_NAME, DATE_ENTERED, INSTRUMENT_SOURCE) '
                            'VALUES (:coverage, :fname, SYSDATE, :instrument)',
                            values)
            con.commit()
            logger.info('Processing successful for DOY %s, year %s ', doy, year)
        except cx_Oracle.DatabaseError as e:
            logger.error('Database entry for DOY %s, year %s has failed: '
                         '%sExiting processing and cleaning that archive file.',
                         doy, year, e)
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

    # remove the files downloaded to the temporary directory
    logger.info('Removing downloaded files')
    System.empty_directory(dloaddir)

    return SUCCESS


def main():
    """ Main routine which grabs the command-line arguments, determines
        which years/days of data need to be processed, then processes the user-
        specified dates of EP/TOMS data.

        Returns:
            SUCCESS/ERROR on success or failure of LADS data retrieval and
            processing
    """
    # get toms config information
    c = AuxConfig()
    cfg = c.get_config('toms')
    start_date = datetime.strptime(cfg.get('toms_start_date'), '%Y-%m-%d')

    # gather the command line arguments
    description = ('Downloads TOMS data and extracts the required parameters '
                   'for Surface Reflectance processing. The parameters are '
                   'then archived for later use.')
    try:
        args = System.get_command_line_arguments(description, start_date, True)
    except Exception as e:
        print(('Argument error: ', e))
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

    # alert user if the overwrite flag is set
    if args.overwrite:
        logger.info('overwrite flag is set: Any existing archive data will '
                    'be overwritten by most current online archive data.')

    # get the command line arguments
    logger.info('Processing EP/TOMS data for %s - %s', args.start_year,
                args.end_year)
    for year in range(args.start_year, args.end_year+1):
        logger.info('Processing year: %s', year)
        status = getTomsData(cfg, year, args.overwrite)
        if status == ERROR:
            logger.warn('Problems occurred while processing EP/TOMS'
                        ' data for year %s.  Processing will continue.', year)

    logger.info('EP/TOMS processing complete.')
    return SUCCESS

if __name__ == "__main__":
    sys.exit(main())
