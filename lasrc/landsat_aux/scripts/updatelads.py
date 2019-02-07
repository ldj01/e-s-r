#!/usr/bin/env python3

import sys
import os
import json
import logging
from datetime import datetime, timedelta

import cx_Oracle
from espa import System
from espa import Date
from espa import HttpSession
from espa import AuxConfig
from api_interface import api_connect

# Global variables
logger = logging.getLogger(__name__)
ERROR = 1
SUCCESS = 0
TERRA_CMA = 0 # Index used for TERRA_CMA files
TERRA_CMG = 1 # Index used for TERRA_CMG files
AQUA_CMA = 2  # Index used for AQUA_CMA files
AQUA_CMG = 3  # Index used for AQUA_CMG files

def downloadLads(year, doy, destination, cfg, session):
    """ Downloads the Aqua/Terra CMG and CMA files for the specified year and
        DOY from the LAADS https interface to the desired destination.
        If the destination directory does not exist, then it is
        made before downloading.  Existing files in the download directory are
        removed/cleaned.

        Parameters:
          year: Year of data to download
          doy: Day of year of data to download
          destination: Name of the directory on the local system to download the
                       LAADS files
          cfg: Config information for filename and archive formats
          session: request session for the desired website

        Returns:
          List of downloaded files in the order [TERRA_CMA, TERRA_CMG, AQUA_CMA,
          AQUA_CMG].  Any missing files are set to None.
    """
    # make sure the download directory exists (and is cleaned up) or create
    # it recursively
    if not os.path.exists(destination):
        logger.info('%s does not exist... creating', destination)
        System.create_directory(destination)
    else:
        logger.info('Cleaning download directory: %s', destination)
        System.empty_directory(destination)

    # obtain the list of URL(s) for our particular date.  this includes the
    # locations for the Aqua and Terra CMG/CMA files.
    urlList = [None] * 4 # create None filled URL list of size 4
    urlList[TERRA_CMA] = cfg.get('terra_cma_format').format(year, doy)
    urlList[TERRA_CMG] = cfg.get('terra_cmg_format').format(year, doy)
    urlList[AQUA_CMA] = cfg.get('aqua_cma_format').format(year, doy)
    urlList[AQUA_CMG] = cfg.get('aqua_cmg_format').format(year, doy)

    # download the data for the specified year from the list of URLs.
    logger.info('Downloading data for %s/%s to %s', year, doy, destination)

    # create a list to keep track of which URLs contain a file
    found_files = []
    for url in urlList:
        try:
            # get the json file associated with this URL.  this is required to
            # determine the filename we want to download
            json_input = json.loads(session.get_lines_from_url(url + '.json')[0]
                                    .decode('utf-8'))
            # there may be more than one file, so sort the file list and only
            # grab the most recent one (ascending sort)
            json_input.sort(key=lambda k: k['last-modified'])
            fname = json_input[-1]['name']
        except Exception:
            logger.warn('No files were found in %s. Continue processing.', url)
            continue

        # download the file
        logger.info('Retrieving %s to %s', os.path.join(url, fname),destination)
        status = session.http_transfer_file(os.path.join(url, fname),
                                            os.path.join(destination, fname))
        if not status:
            found_files.append(None)
            logger.warn('Unsucessful download of %s to %s. Continue processing',
                        (url + fname), os.path.join(destination, fname))
        else:
            found_files.append(os.path.join(destination, fname))
    # end url loop

    return found_files

def getLadsData(cfg, year, session, overwrite):
    """ Downloads the daily MODIS Aqua/Terra CMG and CMA data files for the
        desired year, then combines those files into one daily product
        containing the various ozone, water vapor, temperature, etc. SDSs.

        Parameters:
          cfg: Config information for filename and archive formats
          year: Year of LAADS data to be downloaded and processed
          session: Request session for the desired website
          overwrite: Boolean to determin whether existing archive files are
                     overwritten or skipped

        Returns:
            SUCCESS/ERROR on successful or failed processing
    """
    # gather directory and filename format info
    ancdir = cfg.get('lads_base_archive').format(year)
    dloaddir = cfg.get('lads_temp_directory').format(year)
    outputDir = cfg.get('lads_archive_format').format(ancdir, year)
    archiveFileFormat = cfg.get('lads_filename_format')

    # connect to the database
    try:
        assert os.getenv('IAS_DB_COM') is not None
        con = cx_Oracle.connect(os.getenv('IAS_DB_COM'))
        dbh = con.cursor()
        logger.info('Connected to database successfully')
    except AssertionError as ae:
        logger.error('IAS_DB_COM is not defined')
        return ERROR
    except cx_Oracle.DatabaseError as e:
        logger.error('Could not connect to database: %s', e)
        return ERROR

    # create the output directory if it doesn't exist.
    System.create_directory(outputDir)

    # generate a list of files in the current archive
    archive = System.directory_file_list(outputDir)

    # if the specified year is the current year, only process up through
    # today otherwise process through all the days in the year
    day_of_year = Date.get_days_elapsed_in_year(year)

    # loop through each day in the year and process the LAADS data
    for doy in range(1, day_of_year + 1):
        fname = archiveFileFormat.format(year, doy)
        # skip files that do not need to be overwritten
        if not overwrite and os.path.join(outputDir, fname) in archive:
            continue

        # download the daily LAADS files for the specified year and DOY
        found_files = downloadLads(year, doy, dloaddir, cfg, session)
        logger.debug("Found files: %s", str(found_files))

        # make sure at least one of the Aqua or Terra CMG files is present
        if found_files[TERRA_CMG] is None and found_files[AQUA_CMG] is None:
            logger.warning('No Aqua or Terra LAADS CMG data available for DOY '
                           '%s year %s. Skipping this date.', doy, year)
            continue

        # make sure at least one of the Aqua or Terra CMA files is present
        if found_files[TERRA_CMA] is None and found_files[AQUA_CMA] is None:
            logger.warning('No Aqua or Terra LAADS CMA data available for DOY '
                           '%s year %s. Skipping this date.', doy, year)
            continue

        # prepare cmd line arguments
        terra_cma, terra_cmg, aqua_cma, aqua_cmg = '', '', '', ''
        if found_files[TERRA_CMA] is not None:
            terra_cma = '--terra_cma {}'.format(found_files[TERRA_CMA])
        if found_files[TERRA_CMG] is not None:
            terra_cmg = '--terra_cmg {}'.format(found_files[TERRA_CMG])
        if found_files[AQUA_CMA] is not None:
            aqua_cma = '--aqua_cma {}'.format(found_files[AQUA_CMA])
        if found_files[AQUA_CMG] is not None:
            aqua_cmg = '--aqua_cmg {}'.format(found_files[AQUA_CMG])

        # generate the command-line string combining the CMG and CMA products
        cmdstr = ('combine_l8_aux_data {} {} {} {} --output_dir {} --verbose'
                  .format(terra_cmg, terra_cma, aqua_cmg, aqua_cma, outputDir))
        status = System.execute_cmd(cmdstr)
        if status != 0:
            logger.error('Error running combine_l8_aux_data for DOY %s, year '
                         '%s', year, doy)
            return ERROR

        # remove the files downloaded to the temporary directory
        logger.info('Removing downloaded files from %s', dloaddir)
        System.empty_directory(dloaddir)

        # build the db entry for this file and insert/update the record
        values = {}
        values['coverage'] = datetime(year, 1, 1) + timedelta(days=doy-1)
        values['fname'] = fname
        values['spacecraft'] = 'BOTH'
        if found_files[TERRA_CMA] is None and found_files[TERRA_CMG] is None:
            values['spacecraft'] = 'AQUA'
        if found_files[AQUA_CMA] is None and found_files[AQUA_CMG] is None:
            values['spacecraft'] = 'TERRA'

        try:
            existing = overwrite and os.path.join(outputDir, fname) in archive
            if existing:
                # if we are overwriting and the file is already in the archive,
                # then a db record already exists and must be updated
                logger.info('Updating existing database archive entry')
                dbh.execute('UPDATE MODIS_ATMOS SET FILE_NAME = :fname, '
                            'DATE_ENTERED = SYSDATE, SPACECRAFT_SOURCE = '
                            ':spacecraft WHERE EFFECTIVE_DATE = :coverage',
                            values)
                if dbh.rowcount != 1:
                    logger.error('No rows updated')
                    raise cx_Oracle.DatabaseError
            else:
                logger.info('Creating new database archive entry')
                dbh.execute('INSERT INTO MODIS_ATMOS (EFFECTIVE_DATE, '
                            'FILE_NAME, DATE_ENTERED, SPACECRAFT_SOURCE) '
                            'VALUES (:coverage, :fname, SYSDATE, :spacecraft)',
                            values)
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

    return SUCCESS

def main():
    """ Main routine which grabs the command-line arguments, determines which
        dates and times of data need to be processed, then processes the LADS
        data that falls between the specified dates.

        Returns:
            SUCCESS/ERROR on success or failure of LADS data retrieval and
            processing
    """
    # get lads config information
    c = AuxConfig()
    cfg = c.get_config('lads')
    start_date = datetime.strptime(cfg.get('lads_start_date'), '%Y-%m-%d')

    # Gather command line arguments
    description = ('Downloads LADS data and extracts the required parameters '
                   'for Surface Reflectance processing.  The parameters are '
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

    # Alert user if the overwrite flag is set
    if args.overwrite:
        logger.info('overwrite flag is set: Any existing archive data will '
                    'be overwritten by most current online archive data.')

    # Attempt to get the token from environment
    token = os.getenv('NASA_EARTHDATA_TOKEN')
    # Or from the ESPA API
    if token is None:
        # ESPA Processing Environment
        # Read ~/.usgs/espa/processing.conf to get the URL for the ESPA API.
        # Connect to the ESPA API and get the application token for downloading
        # the LAADS data from the internal database.
        proc_cfg = AuxConfig()
        proc_cfg = proc_cfg.get_config('processing', False)
        server = None
        rpcurl = proc_cfg.get('espa_api')
        if rpcurl is not None:
            server = api_connect(rpcurl)
        if server:
            token = server.get_configuration('aux.downloads.laads.token')

    if token is None:
        logger.error('Application token is None. This needs to be a valid '
                     'token provided for accessing the LAADS data. ')
        return ERROR

    # Create HttpSession object for downloading files and log in to the site
    session = HttpSession()
    session.login(cfg.get('nasa_server_url'), token)

    # Iterate through the given date range and process data
    logger.info('Processing LADS data for %s - %s', args.start_year,
                args.end_year)
    for year in range(args.start_year, args.end_year+1):
        logger.info('Processing year: %s', year)
        status = getLadsData(cfg, year, session, args.overwrite)
        if status == ERROR:
            logger.warn('Problems occurred while processing LADS data for '
                        'year %s. Processing will continue.', year)

    logger.info('LADS processing complete.')
    return SUCCESS


if __name__ == "__main__":
    sys.exit(main())
