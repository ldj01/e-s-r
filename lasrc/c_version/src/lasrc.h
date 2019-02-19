#ifndef _LASRC_H_
#define _LASRC_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"
#include "input.h"
#include "output.h"
#include "lut_subr.h"
#include "espa_metadata.h"
#include "espa_geoloc.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "error_handler.h"

#define ESPA_EPSILON 0.00001

void usage ();

/* data retrieval functions*/
double get_scale_refl();    /* scale for reflective bands */
double get_offset_refl();   /* add offset for reflective bands */
double get_scale_therm();  /* scale for thermal bands */
double get_offset_therm();  /* add offset for thermal bands */
double get_mult_refl();     /* output reflective scale factor */
double get_mult_therm();    /* output thermal scale factor */
int    get_num_threads();   /* number of threads */

int get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **xml_infile,    /* O: address of input XML file */
    char **aux_infile,    /* O: address of input auxiliary file containing
                                water vapor and ozone */
    bool *process_sr,     /* O: process the surface reflectance products */
    bool *write_toa,      /* O: write intermediate TOA products flag */
    bool *verbose         /* O: verbose flag */
);

int compute_toa_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 *sza,         /* I: scaled per-pixel solar zenith angles (degrees),
                              nlines x nsamps */
    float **sband       /* O: output TOA reflectance and brightness temp
                              values (scaled) */
);

int compute_sr_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    char *xml_infile,   /* I: input XML filename */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float pixsize,      /* I: pixel size for the reflectance bands */
    float **sband,      /* I/O: input TOA and output surface reflectance */
    float xts,          /* I: solar zenith angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm         /* I: auxiliary filename for ozone and water vapor */
);

bool is_cloud
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
);

bool is_cloud_or_shadow
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
);

bool is_shadow
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
);

#endif
