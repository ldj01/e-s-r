#!/usr/bin/env python

import sys
import os
import fnmatch
import time
import subprocess
from optparse import OptionParser
import logging

ERROR = 1
SUCCESS = 0

URL = 'http://losrlost02.cr.usgs.gov:7480'


def download_file(url, destination):
    # get the logger
    logger = logging.getLogger(__name__)

    # Download the file.
    # If there'ss a problem with the connection, then retry up to 5 times.
    # Don't use the verbose version which fills the log files with the download
    # percentages.
    msg = 'Retrieving {} to {}'.format(url, destination)
    logger.info(msg)
    cmd = ('wget --tries=5 --no-verbose {}'.format(url))
    retval = subprocess.call(cmd, shell=True, cwd=destination)

    # Make sure the wget was successful, or retry up to 5 more times and
    # sleep in between.
    if retval:
        retry_count = 1
        while (retry_count <= 5 and retval):
            time.sleep(60)
            logger.info('Retry {} of wget for {}'.format(retry_count, url))
            retval = subprocess.call(cmd, shell=True, cwd=destination)
            retry_count += 1

        if retval:
            logger.warn('Unsuccessful download of {} (retried 5 times)'
                        .format(url))

    return SUCCESS


def getLdcmlutData(auxdir):
    # get the logger
    logger = logging.getLogger(__name__)

    # determine the directory for the output auxiliary data files to be
    # processed.  create the directory if it doesn't exist.
    outputDir = '{}/LDCMLUT'.format(auxdir)
    if not os.path.exists(outputDir):
        msg = '{} does not exist... creating'.format(outputDir)
        logger.info(msg)
        os.makedirs(outputDir, 0777)

    file_list = ['AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII',
                 'ANGLE_NEW.hdf',
                 'RES_LUT_V3.0-URBANCLEAN-V2.0.hdf',
                 'TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII',
                 'gascoef-ldcm.ASC',
                 'l8geom.hdf',
                 'tauray-ldcm.ASC']

    for file in file_list:
        # Skip files that already exist.
        skip_file = False
        for myfile in os.listdir(outputDir):
            if fnmatch.fnmatch (myfile, file):
                msg = '{} already exists. Skip.'.format(file)
                logger.info(msg)
                skip_file = True
                break

        if skip_file:
            continue

        # Download the file.
        url = '{}/LDCMLUT/{}'.format(URL, file)
        if download_file(url, outputDir) == ERROR:
            return ERROR

    # Get the LADS files.
    outputDir = '{}/LADS/2018'.format(auxdir)
    if not os.path.exists(outputDir):
        msg = '{} does not exist... creating'.format(outputDir)
        logger.info(msg)
        os.makedirs(outputDir, 0777)

    file_list = ['L8ANC2018001.hdf_fused']

    for file in file_list:
        # Skip files that already exist.
        skip_file = False
        for myfile in os.listdir(outputDir):
            if fnmatch.fnmatch (myfile, file):
                msg = '{} already exists. Skip.'.format(file)
                logger.info(msg)
                skip_file = True
                break

        if skip_file:
            continue

        # Download the file.
        url = '{}/LADS/{}'.format(URL, file)
        if download_file(url, outputDir) == ERROR:
            return ERROR

    # Get the files in the main L8 dir.
    file_list = ['CMGDEM.hdf',
                 'ratiomapndwiexp.hdf']

    for file in file_list:
        # Skip files that already exist.
        skip_file = False
        for myfile in os.listdir(auxdir):
            if fnmatch.fnmatch (myfile, file):
                msg = '{} already exists. Skip.'.format(file)
                logger.info(msg)
                skip_file = True
                break

        if skip_file:
            continue

        # Download the file.
        url = '{}/L8STUFF/{}'.format(URL, file)
        if download_file(url, auxdir) == ERROR:
            return ERROR

    return SUCCESS


def main ():
    logger = logging.getLogger(__name__)  # Get logger for the module.

    # Determine the auxiliary directory to store the data.
    auxdir = os.environ.get('L8_AUX_DIR')
    if auxdir is None:
        msg = 'L8_AUX_DIR environment variable not set... exiting'
        logger.error(msg)
        return ERROR

    status = getLdcmlutData(auxdir)
    if status == ERROR:
        msg = ('Problems occurred while processing LDCMLUT data.')
        logger.error(msg)
        return ERROR

    msg = 'LDCMLUT processing complete.'
    logger.info(msg)
    return SUCCESS


if __name__ == "__main__":
    # setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)
    sys.exit (main())
