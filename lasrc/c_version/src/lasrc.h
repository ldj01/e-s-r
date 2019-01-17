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

/* Defines */
#define ESPA_EPSILON 0.00001

/* Prototypes */
void usage ();

/* data retreval functions*/
double get_scale_refl();    /* scale for reflective bands */
double get_offset_refl();   /* add offset for reflective bands */
double get_scale_therm();  /* scale for thermal bands */
double get_offset_therm();  /* add offset for thermal bands */
double get_mult_refl();     /* scale_refl inverse */
double get_mult_therm();    /* scale_therm inverse */
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

void usage ();

bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is rightmost bit) */
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
    float **sband,      /* O: output TOA reflectance and brightness temp
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
    float **sband,      /* I/O: input TOA and output surface reflectance */
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

int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    float *eps,         /* O: angstrom coefficient */
    int *iaots,         /* O: index for AOTs */
    float *xtv,         /* O: observation zenith angle (deg) */
    float *xmuv,        /* O: cosine of observation zenith angle */
    float *xfi,         /* O: azimuthal difference between sun and
                              observation (deg) */
    float *cosxfi,      /* O: cosine of azimuthal difference */
    float *raot550nm,   /* O: nearest value of AOT */
    float *pres,        /* O: surface pressure */
    float *uoz,         /* O: total column ozone */
    float *uwv,         /* O: total column water vapor (precipital water
                              vapor) */
    float *xtsstep,     /* O: solar zenith step value */
    float *xtsmin,      /* O: minimum solar zenith value */
    float *xtvstep,     /* O: observation step value */
    float *xtvmin,      /* O: minimum observation value */
    float *tsmax,       /* O: maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,       /* O: minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22],      /* O: sun angle table */
    float *ttv,         /* O: view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int32 indts[22],    /* O: index for the sun angle table */
    float *rolutt,      /* O: intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,      /* O: transmission table
                      [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUN_ANGLE_VALS] */
    float *sphalbt,     /* O: spherical albedo table
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,     /* O: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *nbfic,       /* O: communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,        /* O: number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: integer band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: integer band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: integer band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
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
    float band4_pix,     /* I: Band 4 reflectance for current pixel */
    float band5_pix      /* I: Band 5 reflectance for current pixel */
);

bool find_closest_non_fill
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-fill pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-fill pix in aerosol window */
);

bool find_closest_non_cloud_shadow_water
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    float **sband,     /* I: input surface reflectance, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
);

bool find_closest_non_water
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    float **sband,     /* I: input surface reflectance */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
);

void mask_aero_window
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    float **sband,     /* I: input surface reflectance */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    bool *quick_qa     /* O: quick QA for the current aerosol window,
                             AERO_WINDOW x AERO_WINDOW
                             (true=not clear, false=clear) */
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
    if (((l1_qa_pix >> ESPA_L1_DESIGNATED_FILL_BIT) & ESPA_L1_SINGLE_BIT) == 1)
        return true;
    else
        return false;
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
