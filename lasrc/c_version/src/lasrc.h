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
    uint16 **sband,     /* O: output TOA reflectance and brightness temp
                              values (scaled) */
    uint16 *radsat      /* O: radiometric saturation QA band, nlines x nsamps;
                              array should be all zeros on input to this
                              routine*/
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
    uint16 **sband,     /* I/O: input TOA and output surface reflectance */
    int16 *sza,         /* I: per-pixel solar zenith angles, nlines x nsamps */
    int16 *saa,         /* I: per-pixel solar azimuth angles, nlines x nsamps */
    int16 *vza,         /* I: per-pixel view zenith angles, nlines x nsamps */
    int16 *vaa,         /* I: per-pixel view azimuth angles, nlines x nsamps */
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

bool is_water
(
    uint16 band4_pix,     /* I: Band 4 reflectance for current pixel */
    uint16 band5_pix      /* I: Band 5 reflectance for current pixel */
);


/* Defines for the Level-1 BQA band */
/* Define the constants used for shifting bits and ANDing with the bits to
   get to the desire quality bits */
#define ESPA_L1_SINGLE_BIT 0x01             /* 00000001 */
#define ESPA_L1_DOUBLE_BIT 0x03             /* 00000011 */
#define ESPA_L1_DESIGNATED_FILL_BIT 0       /* one bit */
#define ESPA_L1_TERRAIN_OCCLUSION_BIT 1     /* one bit (L8/OLI) */
#define ESPA_L1_RAD_SATURATION_BIT 2        /* two bits */
#define ESPA_L1_CLOUD_BIT 4                 /* one bit */
#define ESPA_L1_CLOUD_CONF_BIT 5            /* two bits */
#define ESPA_L1_CLOUD_SHADOW_CONF_BIT 7     /* two bits */
#define ESPA_L1_SNOW_ICE_CONF_BIT 9         /* two bits */
#define ESPA_L1_CIRRUS_CONF_BIT 11          /* two bits (L8/OLI) */


/******************************************************************************
MODULE:  btest

PURPOSE:  Tests to see if bit n is set in the byte_val variable.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      bit n is not set in byte_val
true       bit n is set in byte_val

NOTES:
******************************************************************************/
static inline bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is least significant
                            bit) */
)
{
    /* Take 2 ** n, then AND that result with the byte value */
    return (byte_val & (1 << n));
}

/******************************************************************************
MODULE:  level1_qa_is_fill

PURPOSE: Determines if the current Level-1 QA pixel is fill

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is fill
false           Pixel is not fill

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool level1_qa_is_fill
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    return ((l1_qa_pix >> ESPA_L1_DESIGNATED_FILL_BIT) & ESPA_L1_SINGLE_BIT)
        == 1;
}

/******************************************************************************
MODULE:  level1_qa_cloud_confidence

PURPOSE: Returns the cloud confidence value (0-3) for the current Level-1 QA
pixel.

RETURN VALUE:
Type = uint8_t
Value           Description
-----           -----------
0               Cloud confidence bits are 00
1               Cloud confidence bits are 01
2               Cloud confidence bits are 10
3               Cloud confidence bits are 11

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline uint8_t level1_qa_cloud_confidence
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    return ((l1_qa_pix >> ESPA_L1_CLOUD_CONF_BIT) & ESPA_L1_DOUBLE_BIT);
}

/******************************************************************************
MODULE:  level1_qa_cloud_shadow_confidence

PURPOSE: Returns the cloud shadow value (0-3) for the current Level-1 QA
pixel.

RETURN VALUE:
Type = uint8_t
Value           Description
-----           -----------
0               Cloud shadow bits are 00
1               Cloud shadow bits are 01
2               Cloud shadow bits are 10
3               Cloud shadow bits are 11

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline uint8_t level1_qa_cloud_shadow_confidence
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    return ((l1_qa_pix >> ESPA_L1_CLOUD_SHADOW_CONF_BIT) & ESPA_L1_DOUBLE_BIT);
}

/******************************************************************************
MODULE:  level1_qa_cirrus_confidence

PURPOSE: Returns the cirrus confidence value (0-3) for the current Level-1 QA
pixel.

RETURN VALUE:
Type = uint8_t
Value           Description
-----           -----------
0               Cirrus confidence bits are 00
1               Cirrus confidence bits are 01
2               Cirrus confidence bits are 10
3               Cirrus confidence bits are 11

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline uint8_t level1_qa_cirrus_confidence
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    return ((l1_qa_pix >> ESPA_L1_CIRRUS_CONF_BIT) & ESPA_L1_DOUBLE_BIT);
}

#endif
