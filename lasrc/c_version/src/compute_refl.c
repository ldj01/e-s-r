#include "lasrc.h"
#include "time.h"
#include "aero_interp.h"
#include "poly_coeff.h"
#include "read_level1_qa.h"

#define DEM_FILL -9999  /* elevation fill value */

/******************************************************************************
MODULE:  compute_toa_refl

PURPOSE:  Computes the TOA reflectance and TOA brightness temps for all
the bands except the pan band. Uses a per-pixel solar zenith angle for the
TOA corrections. Also determines radiometric saturation for each band, as
available.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. These TOA and BT algorithms match those as published by the USGS Landsat
     team in http://landsat.usgs.gov/Landsat8_Using_Product.php
******************************************************************************/
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
                              values (unscaled) */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_toa_refl";   /* function name */
    int i;               /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int iband;           /* current band */
    float refl_mult;     /* reflectance multiplier for bands 1-9 */
    float refl_add;      /* reflectance additive for bands 1-9 */
    float xcals;         /* radiance multiplier for bands 10 and 11 */
    float xcalo;         /* radiance additive for bands 10 and 11 */
    float k1b;           /* K1 temperature constant */
    float k2b;           /* K2 temperature constant */
    uint16 *uband = NULL;  /* array for input image data for a single band,
                              nlines x nsamps */
    time_t mytime;       /* time variable */
    double angband_scale = 1,  /* solar angle band scale factor */
           angband_offset = 0; /* solar angle band offset */


    /* Start the processing */
    mytime = time(NULL);
    printf ("Start TOA reflectance corrections: %s", ctime(&mytime));

    /* Allocate memory for band data */
    uband = malloc (nlines*nsamps*sizeof(uint16));
    if (uband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for uband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the solar angle band scale and offset factors. */
    for (iband = 0; iband < xml_metadata->nbands; iband++)
    {
        Espa_band_meta_t *b = &xml_metadata->band[iband];
        if (strcmp(b->name, "solar_zenith_band4") == 0)
        {
            if (b->scale_factor != ESPA_FLOAT_META_FILL)
                angband_scale = b->scale_factor;
            if (b->add_offset != ESPA_FLOAT_META_FILL)
                angband_offset = b->add_offset;
            break;
        }
    }
    if (iband == xml_metadata->nbands)
    {
        error_handler(true, FUNC_NAME,
                      "Error: Unable to locate solar angle band in metadata.");
        return ERROR;
    }

    /* Loop through all the bands (except the pan band) and compute the TOA
       reflectance and TOA brightness temp */
    for (ib = DN_BAND1; ib <= DN_BAND11; ib++)
    {
        float *band_ptr;  /* convenience pointer */

        /* Read the current band and calibrate bands 1-7 to obtain TOA
           reflectance. Bands are corrected for the sun angle. */
        if (ib <= DN_BAND7)
        {
            iband = ib;
            band_ptr = sband[ib];

            if (get_input_refl_lines (input, iband, 0, nlines, uband) !=
                SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get TOA reflectance coefficients for this reflectance band from
               XML file */
            refl_mult = input->meta.gain[iband];
            refl_add = input->meta.bias[iband];

#ifdef _OPENMP
            #pragma omp parallel for private (i) \
            num_threads(get_num_threads())
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                float xmus;   /* cosine of solar zenith angle (per-pixel) */
                float rotoa;  /* top of atmosphere reflectance */

                /* If this pixel is fill, continue with the next pixel. */
                if (level1_qa_is_fill(qaband[i]))
                {
                    band_ptr[i] = FILL_VALUE;
                    continue;
                }

                /* Compute the TOA reflectance based on the per-pixel
                   sun angle (need to unscale). Scale the TOA value for
                   output. */
                xmus = cos((sza[i]*angband_scale + angband_offset)*DEG2RAD);
                rotoa = (uband[i] * refl_mult) + refl_add;
                rotoa /= xmus;

                /* Save the unscaled TOA reflectance value, but make
                   sure it falls within the defined valid range of
                   unscaled values. */
                if (rotoa < MIN_VALID_REFL)
                    band_ptr[i] = MIN_VALID_REFL;
                else if (rotoa > MAX_VALID_REFL)
                    band_ptr[i] = MAX_VALID_REFL;
                else
                    band_ptr[i] = rotoa;
            }  /* pixel loop */
        }  /* end if band <= band 7 */

        /* Read the current band and calibrate thermal bands.  Not available
           for OLI-only scenes. */
        else if ((ib == DN_BAND10 || ib == DN_BAND11) &&
                 strcmp(instrument, "OLI"))
        {
            int thermal_band_index;  /* 0 = band 10, 1 = band 11 */

            if (ib == DN_BAND10)
            {
                band_ptr = sband[SR_BAND10];
                thermal_band_index = 0;
            }
            else
            {
                band_ptr = sband[SR_BAND11];
                thermal_band_index = 1;
            }

            if (get_input_th_lines(input, thermal_band_index, 0, nlines,
                                   uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+3);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[thermal_band_index];
            xcalo = input->meta.bias_th[thermal_band_index];
            k1b = input->meta.k1_const[thermal_band_index];
            k2b = input->meta.k2_const[thermal_band_index];

            /* Compute brightness temp for band.  Make sure it falls
               within the min/max range for the thermal bands. */

#ifdef _OPENMP
            #pragma omp parallel for private (i) \
            num_threads(get_num_threads())
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                float tmpf;   /* temporary floating point value */

                /* If this pixel is fill, continue with then next pixel. */
                if (level1_qa_is_fill (qaband[i]))
                {
                    band_ptr[i] = FILL_VALUE;
                    continue;
                }

                /* Compute the TOA spectral radiance */
                tmpf = xcals * uband[i] + xcalo;

                /* Compute TOA brightness temp (K). The data
                   will be scaled for output in convert_output() */
                band_ptr[i] = k2b/log(k1b/tmpf + 1.0);
            } /* pixel loop */
        }  /* end if band 10 or 11 */
    }  /* end for ib */

    /* The input data has been read and calibrated. The memory can be freed. */
    free (uband);

    /* Successful completion */
    mytime = time(NULL);
    printf ("End of TOA reflectance computations: %s", ctime(&mytime));
    return (SUCCESS);
}


/******************************************************************************
MODULE:  init_sr_refl

PURPOSE:  Initialization for the atmospheric corrections.  Initialization for
look up tables, auxiliary data, mapping, and geolocation information is used
for the surface reflectance correction.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error initializing the atmospheric parameters
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. The view angle is set to 0.0 and this never changes.
2. The DEM is used to calculate the surface pressure.
******************************************************************************/
static int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid (CMG) DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    float *eps,         /* O: angstrom coefficient */
    float *xtv,         /* O: observation zenith angle (deg) */
    float *xmuv,        /* O: cosine of observation zenith angle */
    float *xfi,         /* O: azimuthal difference between sun and
                              observation (deg) */
    float *cosxfi,      /* O: cosine of azimuthal difference */
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
    int16 *dem,         /* O: climate modeling grid (CMG) DEM data array
                              [DEM_NBLAT x DEM_NBLON] */
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
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "init_sr_refl";       /* function name */
    int retval;          /* return status */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int cmg_pix;         /* pixel location in the CMG array for [lcmg][scmg] */
    int dem_pix;         /* pixel location in the DEM array for [lcmg][scmg] */
    float xcmg, ycmg;    /* x/y location for CMG */

    /* Vars for forward/inverse mapping space */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
    float center_lat, center_lon; /* lat/long for scene center */

    /* Initialize the look up tables */
    *eps = 1.0;
    *xtv = 0.0;
    *xmuv = cos (*xtv * DEG2RAD);
    *xfi = 0.0;
    *cosxfi = cos (*xfi * DEG2RAD);
    *xtsmin = 0;
    *xtsstep = 4.0;
    *xtvmin = 2.84090;
    *xtvstep = 6.52107 - *xtvmin;
    retval = readluts (tsmax, tsmin, ttv, tts, nbfic, nbfi, indts, rolutt,
        transt, sphalbt, normext, *xtsstep, *xtsmin, anglehdf,
        intrefnm, transmnm, spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("The LUTs for urban clean case v2.0 have been read.  We can "
        "now perform atmospheric correction.\n");

    /* Read the auxiliary data files used as input to the reflectance
       calculations */
    retval = read_auxiliary_files (cmgdemnm, rationm, auxnm, dem, andwi, sndwi,
        ratiob1, ratiob2, ratiob7, intratiob1, intratiob2, intratiob7,
        slpratiob1, slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the auxiliary files");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Getting parameters for atmospheric correction */
    /* Update to get the parameter of the scene center */
    *pres = ATMOS_PRES_0;
    *uoz = 0.30;
    *uwv = 0.5;

    /* Use scene center (and center of the pixel) to compute atmospheric
       parameters */
    img.l = (int) (nlines * 0.5) - 0.5;
    img.s = (int) (nsamps * 0.5) + 0.5;
    img.is_fill = false;
    if (!from_space (space, &img, &geo))
    {
        sprintf (errmsg, "Mapping scene center to geolocation coords");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    center_lat = geo.lat * RAD2DEG;
    center_lon = geo.lon * RAD2DEG;
    printf ("Scene center line/sample: %f, %f\n", img.l, img.s);
    printf ("Scene center lat/long: %f, %f\n", center_lat, center_lon);

    /* Use the scene center lat/long to determine the line/sample in the
       CMG-related lookup tables, using the center of the UL pixel.
       Negative latitude values should be the largest line values in the CMG
       grid.  Negative longitude values should be the smallest sample values
       in the CMG grid. */
    ycmg = (89.975 - center_lat) * 20.0;    /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;   /* vs / 0.05 */
    lcmg = (int) roundf (ycmg);
    scmg = (int) roundf (xcmg);

    /* Handle the edges of the lat/long values in the CMG grid */
    if (lcmg < 0)
        lcmg = 0;
    else if (lcmg >= CMG_NBLAT)
        lcmg = CMG_NBLAT - 1;

    if (scmg < 0)
        scmg = 0;
    else if (scmg >= CMG_NBLON)
        scmg = CMG_NBLON - 1;

    cmg_pix = lcmg * CMG_NBLON + scmg;
    if (wv[cmg_pix] != 0)
        *uwv = wv[cmg_pix] / 200.0;
    else
        *uwv = 0.5;

    if (oz[cmg_pix] != 0)
        *uoz = oz[cmg_pix] / 400.0;
    else
        *uoz = 0.3;

    dem_pix = lcmg * DEM_NBLON + scmg;
    if (dem[dem_pix] != DEM_FILL)
        *pres = ATMOS_PRES_0 * exp (-dem[dem_pix] * ONE_DIV_8500);
    else
        *pres = ATMOS_PRES_0;

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  find_closest_non_fill

PURPOSE:  Finds the closest non-fill pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
static bool find_closest_non_fill
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-fill pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-fill pix in aerosol window */
)
{
    int line, samp;          /* looping variables for lines and samples */
    int aero_window;         /* looping variable for the aerosol window */
    uint16 *band_ptr;        /* band data pointer */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        /* Find the starting and ending samples of the window, clamping to
           the image boundaries. */
        int start_samp, end_samp;
        if (center_samp - aero_window > 0)
            start_samp = center_samp - aero_window;
        else
            start_samp = 0;
        if (center_samp + aero_window < nsamps)
            end_samp = center_samp + aero_window;
        else
            end_samp = nsamps - 1;

        /* Check top line of current window. */
        line = center_line - aero_window;
        if (line >= 0)
        {
            band_ptr = &qaband[line*nsamps];
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(band_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
            line++;
            band_ptr += nsamps;
        }
        else
        {
            line = 0;
            band_ptr = qaband;
        }

        /* Check first and last samples of the sides of the current window. */
        for (; line < center_line + aero_window; line++, band_ptr+=nsamps)
        {
            /* If we hit the end of the image, we're done checking the
               sides. */
            if (line >= nlines)
                break;

            samp = center_samp - aero_window;
            if (samp >= 0 && !level1_qa_is_fill(band_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }

            samp = center_samp + aero_window;
            if (samp < nsamps && !level1_qa_is_fill(band_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }
        }

        /* Check the bottom line of the current window. */
        if (line < nlines)
        {
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(band_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
        }
    } /* window loop */

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  find_closest_non_cloud_shadow_water

PURPOSE:  Finds the closest non-cloud, non-shadow, non-water pixel in the
aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
static bool find_closest_non_cloud_shadow_water
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
)
{
    int line, samp;          /* looping variables for lines and samples */
    int line_index;
    int aero_window;         /* looping variable for the aerosol window */
    uint16 *qband_ptr;       /* band data pointer */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        /* Find the starting and ending samples of the window, clamping to
           the image boundaries. */
        int start_samp, end_samp;
        if (center_samp - aero_window > 0)
            start_samp = center_samp - aero_window;
        else
            start_samp = 0;
        if (center_samp + aero_window < nsamps)
            end_samp = center_samp + aero_window;
        else
            end_samp = nsamps - 1;

        /* Check top line of current window. */
        line = center_line - aero_window;
        if (line >= 0)
        {
            line_index = line*nsamps;
            qband_ptr = &qaband[line_index];
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(qband_ptr[samp]) &&
                    !is_cloud_or_shadow(qband_ptr[samp]) &&
                    !level1_qa_is_water(qband_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
            line++;
            qband_ptr += nsamps;
        }
        else
        {
            line = 0;
            qband_ptr = qaband;
        }

        /* Check first and last samples of the sides of the current window. */
        for (; line < center_line + aero_window; line++, qband_ptr+=nsamps)
        {
            /* If we hit the end of the image, we're done checking the
               sides. */
            if (line >= nlines)
                break;

            samp = center_samp - aero_window;
            if (samp >= 0 &&
                !level1_qa_is_fill(qband_ptr[samp]) &&
                !is_cloud_or_shadow(qband_ptr[samp]) &&
                !level1_qa_is_water(qband_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }

            samp = center_samp + aero_window;
            if (samp < nsamps &&
                !level1_qa_is_fill(qband_ptr[samp]) &&
                !is_cloud_or_shadow(qband_ptr[samp]) &&
                !level1_qa_is_water(qband_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }
        }

        /* Check the bottom line of the current window. */
        if (line < nlines)
        {
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(qband_ptr[samp]) &&
                    !is_cloud_or_shadow(qband_ptr[samp]) &&
                    !level1_qa_is_water(qband_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
        }
    } /* window loop */

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  find_closest_non_water

PURPOSE:  Finds the closest non-water pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
static bool find_closest_non_water
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
)
{
    int line, samp;          /* looping variables for lines and samples */
    int line_index;
    int aero_window;         /* looping variable for the aerosol window */
    uint16 *qband_ptr;       /* band data pointer */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        /* Find the starting and ending samples of the window, clamping to
           the image boundaries. */
        int start_samp, end_samp;
        if (center_samp - aero_window > 0)
            start_samp = center_samp - aero_window;
        else
            start_samp = 0;
        if (center_samp + aero_window < nsamps)
            end_samp = center_samp + aero_window;
        else
            end_samp = nsamps - 1;

        /* Check top line of current window. */
        line = center_line - aero_window;
        if (line >= 0)
        {
            line_index = line*nsamps;
            qband_ptr = &qaband[line_index];
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(qband_ptr[samp]) &&
                    !level1_qa_is_water(qband_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
            line++;
            qband_ptr += nsamps;
        }
        else
        {
            line = 0;
            qband_ptr = qaband;
        }

        /* Check first and last samples of the sides of the current window. */
        for (; line < center_line + aero_window; line++, qband_ptr+=nsamps)
        {
            /* If we hit the end of the image, we're done checking the
               sides. */
            if (line >= nlines)
                break;

            samp = center_samp - aero_window;
            if (samp >= 0 &&
                !level1_qa_is_fill(qband_ptr[samp]) &&
                !level1_qa_is_water(qband_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }

            samp = center_samp + aero_window;
            if (samp < nsamps &&
                !level1_qa_is_fill(qband_ptr[samp]) &&
                !level1_qa_is_water(qband_ptr[samp]))
            {
                *nearest_line = line;
                *nearest_samp = samp;
                return true;
            }
        }

        /* Check the bottom line of the current window. */
        if (line < nlines)
        {
            for (samp = start_samp; samp <= end_samp; samp++)
            {
                if (!level1_qa_is_fill(qband_ptr[samp]) &&
                    !level1_qa_is_water(qband_ptr[samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return true;
                }
            }
        }
    } /* window loop */

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  compute_sr_refl

PURPOSE:  Computes the surfance reflectance for all the reflectance bands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Initializes the variables and data arrays from the lookup table and
   auxiliary files.
2. The tauray array was originally read in from a static ASCII file, but it is
   now hardcoded to save time from reading the file each time.  This file was
   generated (like many of the other auxiliary input tables) by running 6S and
   storing the coefficients.
4. Aerosols are retrieved for all non-fill pixels.  If the aerosol fails the
   model residual or NDVI test, then the pixel is flagged as water.  All water
   pixels are run through a water-specific aerosol retrieval.  If the model
   residual fails, then that pixel is marked as failed aerosol retrieval.  Any
   pixel that failed retrieval is then interpolated using an average of the
   clear (valid land pixel aerosols) and water (valid water pixel aerosols).
   Those final aerosol values are used for the surface reflectance corrections.
5. Cloud-based QA information is not processed in this algorithm.
******************************************************************************/
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
    float xts,          /* I: scene center solar zenith angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid (CMG) DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm         /* I: auxiliary filename for ozone and water vapor */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_sr_refl";   /* function name */
    int retval;          /* return status */
    int i;               /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int center_line;     /* aerosol window center line */
    float tgo;           /* other gaseous transmittance (tgog * tgoz) */
    float roatm;         /* intrinsic atmospheric reflectance */
    float ttatmg;        /* total atmospheric transmission */
    float satm;          /* atmosphere spherical albedo */
    float xrorayp;       /* reflectance of the atmosphere due to molecular
                            (Rayleigh) scattering */
    float next;
    float btgo[NSR_BANDS];     /* other gaseous transmittance for bands 1-7 */
    float broatm[NSR_BANDS];   /* atmospheric reflectance for bands 1-7 */
    float bttatmg[NSR_BANDS];  /* ttatmg for bands 1-7 */
    float bsatm[NSR_BANDS];    /* atmosphere spherical albedo for bands 1-7 */
    float median_aerosol; /* median aerosol value for clear pixels */
    uint8 *ipflag = NULL; /* QA flag to assist with aerosol interpolation,
                             nlines x nsamps */
    float *taero = NULL;  /* aerosol values for each pixel, nlines x nsamps */
    float *teps = NULL;   /* angstrom coeff for each pixel, nlines x nsamps */
    float *aerob1 = NULL; /* atmospherically corrected band 1 data
                             (TOA refl), nlines x nsamps */
    float *aerob2 = NULL; /* atmospherically corrected band 2 data
                             (TOA refl), nlines x nsamps */
    float *aerob4 = NULL; /* atmospherically corrected band 4 data
                             (TOA refl), nlines x nsamps */
    float *aerob5 = NULL; /* atmospherically corrected band 5 data
                             (TOA refl), nlines x nsamps */
    float *aerob7 = NULL; /* atmospherically corrected band 7 data
                             (TOA refl), nlines x nsamps */
    uint16 *out_band = NULL;  /* scaled output */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */

    /* Lookup table variables */
    float eps;           /* angstrom coefficient */
    float xtv;           /* observation zenith angle (deg) */
    float xmuv;          /* cosine of observation zenith angle */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float cosxfi;        /* cosine of azimuthal difference */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */
    float *rolutt = NULL;  /* intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt = NULL;  /* transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float *sphalbt = NULL; /* spherical albedo table 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext = NULL; /* aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax = NULL;   /* maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin = NULL;   /* minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi = NULL;    /* number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic = NULL;   /* communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *ttv = NULL;     /* view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22];         /* sun angle table */
    int32 indts[22];       /* index for sun angle table */

    /* Atmospheric correction coefficient variables */
    float tgo_arr[NREFL_BANDS];     /* per-band other gaseous transmittance */
    float xrorayp_arr[NREFL_BANDS]; /* per-band reflectance of the atmosphere
                                       due to molecular (Rayleigh) scattering */
    float roatm_arr[NREFL_BANDS][NAOT_VALS];  /* per band AOT vals for roatm */
    float ttatmg_arr[NREFL_BANDS][NAOT_VALS]; /* per band AOT vals for ttatmg */
    float satm_arr[NREFL_BANDS][NAOT_VALS];   /* per band AOT vals for satm */
    float roatm_coef[NREFL_BANDS][NCOEF];  /* per band poly coeffs for roatm */
    float ttatmg_coef[NREFL_BANDS][NCOEF]; /* per band poly coeffs for ttatmg */
    float satm_coef[NREFL_BANDS][NCOEF];   /* per band poly coeffs for satm */
    float normext_p0a3_arr[NREFL_BANDS];   /* per band normext[iband][0][3] */
    int roatm_iaMax[NREFL_BANDS];          /* ??? */
    int ia;                                /* looping variable for AOTs */
    int iaMaxTemp;                         /* max temp for current AOT level */

    /* Auxiliary file variables */
    int16 *dem = NULL;        /* CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi = NULL;      /* avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi = NULL;      /* standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1 = NULL;    /* mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2 = NULL;    /* mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7 = NULL;    /* mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1 = NULL;   /* intercept band1 ratio,
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob2 = NULL;   /* intercept band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob7 = NULL;   /* intercept band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob1 = NULL;   /* slope band1 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob2 = NULL;   /* slope band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob7 = NULL;   /* slope band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    uint16 *wv = NULL;       /* water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz = NULL;        /* ozone values [CMG_NBLAT x CMG_NBLON] */
    float uoz;          /* total column ozone */
    float uwv;          /* total column water vapor (precipital water vapor) */
    float pres;         /* surface pressure */

    /* Output file info */
    time_t mytime;               /* timing variable */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    char envi_file[STR_SIZE];    /* ENVI filename */
    char *cptr = NULL;           /* pointer to the file extension */

    /* Table constants */
    float aot550nm[NAOT_VALS] =  /* AOT look-up table */
        {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.20,
         1.40, 1.60, 1.80, 2.00, 2.30, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00};
    float tpres[NPRES_VALS] =    /* surface pressure table */
        {1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0};

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities */
    float tauray[NSR_BANDS] =  /* molecular optical thickness coefficients --
        produced by running 6S */
        {0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037,
         0.07984};
    double oztransa[NSR_BANDS] =   /* ozone transmission coeff */
        {-0.00255649, -0.0177861, -0.0969872, -0.0611428, 0.0001, 0.0001,
          0.0001, -0.0834061};
    double wvtransa[NSR_BANDS] =   /* water vapor transmission coeff */
        {2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159, 0.000729136,
         0.00067324, 0.0177533, 0.00279738};
    double wvtransb[NSR_BANDS] =   /* water vapor transmission coeff */
        {0.999742, 0.999742, 0.775024, 0.774482, 0.893085, 0.939669, 0.65094,
         0.759952};
    double ogtransa1[NSR_BANDS] =  /* other gases transmission coeff */
        {4.91586e-20, 4.91586e-20, 4.91586e-20, 1.04801e-05, 1.35216e-05,
         0.0205425, 0.0256526, 0.000214329};
    double ogtransb0[NSR_BANDS] =  /* other gases transmission coeff */
        {0.000197019, 0.000197019, 0.000197019, 0.640215, -0.195998, 0.326577,
         0.243961, 0.396322};
    double ogtransb1[NSR_BANDS] =  /* other gases transmission coeff */
        {9.57011e-16, 9.57011e-16, 9.57011e-16, -0.348785, 0.275239, 0.0117192,
         0.0616101, 0.04728};

#ifdef WRITE_TAERO
    FILE *aero_fptr=NULL;   /* file pointer for aerosol files */
#endif

    /* Start processing */
    mytime = time(NULL);
    printf ("Start surface reflectance corrections: %s", ctime(&mytime));

    /* Allocate memory for the many arrays needed to do the surface reflectance
       computations */
    retval = memory_allocation_sr (nlines, nsamps, &aerob1, &aerob2, &aerob4,
        &aerob5, &aerob7, &ipflag, &taero, &teps, &dem,
        &andwi, &sndwi, &ratiob1, &ratiob2, &ratiob7, &intratiob1, &intratiob2,
        &intratiob7, &slpratiob1, &slpratiob2, &slpratiob7, &wv, &oz, &rolutt,
        &transt, &sphalbt, &normext, &tsmax, &tsmin, &nbfic, &nbfi, &ttv,
        &out_band);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error allocating memory for the data arrays needed "
            "for surface reflectance calculations.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the geolocation space applications */
    if (!get_geoloc_info (xml_metadata, &space_def))
    {
        sprintf (errmsg, "Getting the space definition from the XML file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    space = setup_mapping (&space_def);
    if (space == NULL)
    {
        sprintf (errmsg, "Setting up the geolocation mapping");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the look up tables and atmospheric correction variables
       view zenith initialized to 0.0 (xtv)
       azimuthal difference between sun and obs angle initialize to 0.0 (xfi)
       surface pressure is initialized to the pressure at the center of the
           scene (using the DEM) (pres)
       water vapor is initialized to the value at the center of the scene (uwv)
       ozone is initialized to the value at the center of the scene (uoz) */
    retval = init_sr_refl (nlines, nsamps, input, space, anglehdf, intrefnm,
        transmnm, spheranm, cmgdemnm, rationm, auxnm, &eps, &xtv,
        &xmuv, &xfi, &cosxfi, &pres, &uoz, &uwv, &xtsstep,
        &xtsmin, &xtvstep, &xtvmin, tsmax, tsmin, tts, ttv, indts, rolutt,
        transt, sphalbt, normext, nbfic, nbfi, dem, andwi, sndwi, ratiob1,
        ratiob2, ratiob7, intratiob1, intratiob2, intratiob7, slpratiob1,
        slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error initializing the lookup tables and "
            "atmospheric correction variables.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the reflectance bands and perform atmospheric
       corrections based on climatology */
    mytime = time(NULL);
    printf ("Performing atmospheric corrections for each reflectance "
            "band.  %s", ctime(&mytime));
    for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
    {
        float *source_ptr = sband[ib];/* TOA band pointer */
        float *band_ptr;               /* convenience pointer */
        float rotoa;                   /* top of atmosphere reflectance */
        float roslamb;                 /* lambertian surface reflectance */

        /* Get the parameters for the atmospheric correction */
        /* rotoa is not defined for this call, which is ok, but the
           roslamb value is not valid upon output. Just set it to 0.0 to
           be consistent. */
        rotoa = 0.0;
        eps = 2.5;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
            aot550nm[1], ib, pres, tpres, aot550nm, rolutt, transt, xtsstep,
            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
            nbfi, tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
            ogtransb1, wvtransa, wvtransb, oztransa, rotoa, &roslamb,
            &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Save these band-related parameters for later */
        btgo[ib] = tgo;
        broatm[ib] = roatm;
        bttatmg[ib] = ttatmg;
        bsatm[ib] = satm;

        /* Before completing atmospheric corrections, swap the band pointers
           in order to retain the unscaled TOA reflectance values for later
           use. */
        if (ib == DN_BAND1)
        {
            sband[ib] = aerob1;
            aerob1 = source_ptr;
        }
        else if (ib == DN_BAND2)
        {
            sband[ib] = aerob2;
            aerob2 = source_ptr;
        }
        else if (ib == DN_BAND4)
        {
            sband[ib] = aerob4;
            aerob4 = source_ptr;
        }
        else if (ib == DN_BAND5)
        {
            sband[ib] = aerob5;
            aerob5 = source_ptr;
        }
        else if (ib == DN_BAND7)
        {
            sband[ib] = aerob7;
            aerob7 = source_ptr;
        }
        band_ptr = sband[ib];

        /* Perform atmospheric corrections for bands 1-7 */
#ifdef _OPENMP
        #pragma omp parallel for private (i) \
        num_threads(get_num_threads())
#endif
        for (i = 0; i < nlines*nsamps; i++)
        {
            float roslamb; /* lambertian surface reflectance */

            /* Skip fill pixels, which have already been marked in the
               TOA calculations. */
            if (level1_qa_is_fill(qaband[i]))
                continue;

            /* Apply the atmospheric corrections (ignoring the Rayleigh
               scattering component and water vapor), and store the
               unscaled value for further corrections.  (NOTE: the full
               computations are in atmcorlamb2) */
            roslamb = source_ptr[i] - tgo*roatm;
            roslamb /= tgo*ttatmg + satm*roslamb;
            band_ptr[i] = roslamb;
        }  /* pixel loop */
    }  /* for ib */

    /* Start the retrieval of atmospheric correction parameters for each band */
    mytime = time(NULL);
    printf ("Starting retrieval of atmospheric correction parameters.  %s",
        ctime(&mytime));
    for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
    {
        float rotoa;   /* top of atmosphere reflectance */
        float roslamb; /* lambertian surface reflectance */

        /* Get the parameters for the atmospheric correction */
        /* rotoa is not defined for this call, which is ok, but the
           roslamb value is not valid upon output. Just set it to 0.0 to
           be consistent. */
        normext_p0a3_arr[ib] = normext[ib * NPRES_VALS * NAOT_VALS + 0 + 3];
            /* normext[ib][0][3]; */
        rotoa = 0.0;
        eps = 2.5;
        for (ia = 0; ia < NAOT_VALS; ia++)
        {
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                aot550nm[ia], ib, pres, tpres, aot550nm, rolutt, transt,
                xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
                rotoa, &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next,
                eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2 for band %d.", ib);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Store the AOT-related variables for use in the atmospheric
               corrections */
            roatm_arr[ib][ia] = roatm;
            ttatmg_arr[ib][ia] = ttatmg;
            satm_arr[ib][ia] = satm;
        }

        /* Store the band-related variables for use in the atmospheric
           corrections. tgo and xrorayp are the same for each AOT, so just
           save the last set for this band. */
        tgo_arr[ib] = tgo;
        xrorayp_arr[ib] = xrorayp;
    }

    for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
    {
        iaMaxTemp = 1;

        for (ia = 1; ia < NAOT_VALS; ia++)
        {
            if (ia == NAOT_VALS-1)
                iaMaxTemp = NAOT_VALS-1;

            if ((roatm_arr[ib][ia] - roatm_arr[ib][ia-1]) > ESPA_EPSILON)
                continue;
            else
            {
                iaMaxTemp = ia-1;
                break;
            }
        }

        roatm_iaMax[ib] = iaMaxTemp;
        get_3rd_order_poly_coeff(aot550nm, roatm_arr[ib], iaMaxTemp,
                                 roatm_coef[ib]);

        /* Get the polynomial coefficients for ttatmg */
        get_3rd_order_poly_coeff(aot550nm, ttatmg_arr[ib], NAOT_VALS,
                                 ttatmg_coef[ib]);

        /* Get the polynomial coefficients for satm */
        get_3rd_order_poly_coeff(aot550nm, satm_arr[ib], NAOT_VALS,
                                 satm_coef[ib]);
    }

    /* Start the aerosol inversion */
    mytime = time(NULL);
    printf ("Aerosol Inversion using %d x %d aerosol window  %s",
        AERO_WINDOW, AERO_WINDOW, ctime(&mytime));
#ifdef _OPENMP
    #pragma omp parallel for private (center_line) \
    num_threads(get_num_threads())
#endif
    for (center_line = HALF_AERO_WINDOW; center_line < nlines;
         center_line += AERO_WINDOW)
    {
        int center_samp; /* aerosol window center sample */
        int center_pix = center_line*nsamps + HALF_AERO_WINDOW;
               /* current pixel in 1D arrays of nlines*nsamps for the center
                  of the aerosol window */

        for (center_samp = HALF_AERO_WINDOW; center_samp < nsamps;
             center_samp += AERO_WINDOW, center_pix += AERO_WINDOW)
        {
            Img_coord_float_t img; /* coordinate in line/sample space */
            Geo_coord_t geo;  /* coordinate in lat/long space */
            float lat, lon;   /* pixel lat, long location */
            float xcmg, ycmg; /* x/y location for climate modeling grid (CMG) */
            int lcmg, scmg;   /* line/sample index for the CMG */
            int lcmg1, scmg1; /* line+1/sample+1 index for the CMG */
            float u, v;       /* line/sample index for the CMG */
            float u_x_v;      /* u * v */
            int ratio_pix[4]; /* pixel loc for ratio products [lcmg][scmg],
                                 [lcmg][scmg+1], [lcmg+1][scmg],
                                 [lcmg+1][scmg+1] */
            float slpr[3][4];    /* band ratio slopes at (line,samp),
                                    (line,samp+1), (line+1,samp), and
                                    (line+1,samp+1) for bands 1, 2, and 7 */
            float intr[3][4];    /* band ratio intercepts at (line,samp),
                                    (line,samp+1), (line+1,samp), and
                                    (line+1,samp+1) for bands 1, 2, and 7 */
            int ratio_index;  /* index for ratio arrays */
            float slprb[3];   /* interpolated band ratio slope values for
                                 band ratios 1, 2, 7 */
            float intrb[3];   /* interpolated band ratio intercept values for
                                 band ratios 1, 2, 7 */
            int band_index;   /* index for band arrays */
            float xndwi;      /* calculated NDWI value */
            float ndwi_th1, ndwi_th2;     /* values for NDWI calculations */
            int iband;        /* current band */
            int iband1, iband3;           /* band indices (zero-based) */
            float eps1, eps2, eps3;       /* eps values for three runs */
            float epsmin;     /* eps which minimizes the residual */
            double xa, xb;    /* coefficients for finding the eps that
                                 minimizes the residual*/
            float rotoa;      /* top of atmosphere reflectance */
            float residual;   /* model residual */
            float residual1, residual2, residual3; /* residuals for 3 different
                                                      eps values */
            float raot;       /* AOT reflectance */
            float sraot1, sraot3;         /* raot values for three different
                                             eps values */
            float corf;       /* aerosol impact (higher values represent high
                                 aerosol) */
            float ros4, ros5; /* surface reflectance for bands 4 and 5 */
            float erelc[NSR_BANDS];       /* band ratio variable for
                                             bands 1-7 */
            float troatm[NSR_BANDS];      /* atmospheric reflectance table for
                                             bands 1-7 */
            int iaots;        /* index for AOTs */
            int nearest_line = center_line;  /* line for nearest non-fill/cloud
                                                pixel in the aerosol window */
            int nearest_samp = center_samp;  /* samp for nearest non-fill/cloud
                                                pixel in the aerosol window */
            int nearest_pix = center_pix; /* nearest non-fill/cloud pixel in
                                             the aerosol window */

            /* If this pixel is fill */
            if (level1_qa_is_fill (qaband[nearest_pix]))
            {
                /* Look for other non-fill pixels in the window */
                if (find_closest_non_fill (qaband, nlines, nsamps,
                                           center_line, center_samp,
                                           &nearest_line, &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill pixel for
                       further processing of aerosols. However we will still
                       write to the center of the aerosol window for the
                       current window. */
                    nearest_pix = nearest_line*nsamps + nearest_samp;
                }
                else
                {
                    /* No other non-fill pixels found.  Pixel is already
                       flagged as fill. Move to next aerosol window. */
                    continue;
                }
            }

            /* If this non-fill pixel is water, then look for a pixel which is
               not water.  If none are found then the whole window is fill or
               water.  Flag this pixel as water. */
            if (level1_qa_is_water (qaband[nearest_pix]))
            {
                /* Look for other non-fill/non-water pixels in the window.
                   Start with the center of the window and search outward. */
                if (find_closest_non_water (qaband, nlines, nsamps,
                                            center_line, center_samp,
                                            &nearest_line, &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill/non-water
                       pixel for further processing */
                    nearest_pix = nearest_line*nsamps + nearest_samp;
                }
                else
                {
                    /* Assign generic values for the water pixel */
                    ipflag[center_pix] = (1 << IPFLAG_WATER);
                    taero[center_pix] = DEFAULT_AERO;
                    teps[center_pix] = DEFAULT_EPS;

                    /* Next window */
                    continue;
                }
            }

            /* If this non-fill/non-water pixel is cloud or shadow, then look
               for a pixel which is not cloudy, shadow, water, or fill.  If
               none are found, then just use this pixel. */
            if (is_cloud_or_shadow (qaband[nearest_pix]))
            {
                /* Look for other non-fill/non-water/non-cloud/non-shadow
                   pixels in the window.  Start with the center of the window
                   and search outward. */
                if (find_closest_non_cloud_shadow_water (qaband, nlines,
                    nsamps, center_line, center_samp, &nearest_line,
                    &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill/non-cloud
                       pixel for further processing */
                    nearest_pix = nearest_line*nsamps + nearest_samp;
                }
                else
                {
                    /* Assign generic aerosol values for the cloud or shadow
                       pixel. */
                    if (is_cloud (qaband[nearest_pix]))
                        ipflag[center_pix] = (1 << IPFLAG_CLOUD);
                    else if (is_shadow (qaband[nearest_pix]))
                        ipflag[center_pix] = (1 << IPFLAG_SHADOW);
                    taero[center_pix] = DEFAULT_AERO;
                    teps[center_pix] = DEFAULT_EPS;

                    /* Next window */
                    continue;
                }
            }

            /* Get the lat/long for the center of the current pixel
               (which may not be the center of the aerosol window). */
            img.l = nearest_line - 0.5;
            img.s = nearest_samp + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                    "geolocation coords", nearest_line, nearest_samp);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /* Use that lat/long to determine the line/sample in the
               climate modeling grid (CMG)-related lookup tables, using the
               center of the UL pixel.  Note, we are basically making sure
               the line/sample combination falls within -90, 90 and -180, 180
               global climate data boundaries.  However, the source code
               below uses lcmg+1 and scmg+1, which for some scenes may wrap
               around the dateline or the poles.  Thus we need to wrap the
               CMG data around to the beginning of the array. */
            /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
               pixel for each calculation.  Negative latitude values should
               be the largest line values in the CMG grid.  Negative
               longitude values should be the smallest sample values in the
               CMG grid. */
            /* The line/sample calculation from the x/ycmg values are not
               rounded.  The interpolation of the value using line+1 and
               sample+1 are based on the truncated numbers, therefore
               rounding up is not appropriate. */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) ycmg;
            scmg = (int) xcmg;

            /* Handle the edges of the lat/long values in the CMG grid */
            if (lcmg < 0)
                lcmg = 0;
            else if (lcmg >= CMG_NBLAT)
                lcmg = CMG_NBLAT - 1;

            if (scmg < 0)
                scmg = 0;
            else if (scmg >= CMG_NBLON)
                scmg = CMG_NBLON - 1;

            /* If the current CMG pixel is at the edge of the CMG array, then
               allow the next pixel for interpolation to wrap around the
               array */
            if (scmg >= CMG_NBLON - 1)  /* 180 degrees so wrap around */
                scmg1 = 0;
            else
                scmg1 = scmg + 1;

            if (lcmg >= CMG_NBLAT - 1)  /* -90 degrees, so set the next pixel
                                           to also use -90. */
                lcmg1 = lcmg;
            else
                lcmg1 = lcmg + 1;

            /* Determine the fractional difference between the integer location
               and floating point pixel location to be used for interpolation */
            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            u_x_v = u * v;

            /* Determine the band ratios and slope/intercept */
            ratio_pix[0] = lcmg*RATIO_NBLON + scmg;
            ratio_pix[1] = lcmg*RATIO_NBLON + scmg1;
            ratio_pix[2] = lcmg1*RATIO_NBLON + scmg;
            ratio_pix[3] = lcmg1*RATIO_NBLON + scmg1;

            for (ratio_index = 0; ratio_index < 4; ratio_index++)
            {
                int rindex = ratio_pix[ratio_index]; /* convenience variable */

                if (ratiob2[rindex] > 1000 || ratiob1[rindex] > 1000 ||
                    ratiob2[rindex] < 100 || ratiob1[rindex] < 100)
                                                   /* values scaled by 1000 */
                {
                    slpr[0][ratio_index] = 0;
                    slpr[1][ratio_index] = 0;
                    slpr[2][ratio_index] = 0;
                    intr[0][ratio_index] = 550;
                    intr[1][ratio_index] = 600;
                    intr[2][ratio_index] = 2000;
                }
                else if (sndwi[rindex] < 200)
                {
                    slpr[0][ratio_index] = 0;
                    slpr[1][ratio_index] = 0;
                    slpr[2][ratio_index] = 0;
                    intr[0][ratio_index] = ratiob1[rindex];
                    intr[1][ratio_index] = ratiob2[rindex];
                    intr[2][ratio_index] = ratiob7[rindex];
                }
                else
                {
                    slpr[0][ratio_index] = slpratiob1[rindex];
                    slpr[1][ratio_index] = slpratiob2[rindex];
                    slpr[2][ratio_index] = slpratiob7[rindex];
                    intr[0][ratio_index] = intratiob1[rindex];
                    intr[1][ratio_index] = intratiob2[rindex];
                    intr[2][ratio_index] = intratiob7[rindex];
                }
            }

            /* Interpolate the slope/intercept for each band (1, 2, and 7),
               and unscale */
            for (band_index = 0; band_index < 3; band_index++)
            {
                slprb[band_index] = slpr[band_index][0]
                    + u*(slpr[band_index][2] - slpr[band_index][0])
                    + v*(slpr[band_index][1] - slpr[band_index][0])
                    + u_x_v*(slpr[band_index][0] - slpr[band_index][1]
                             - slpr[band_index][2] + slpr[band_index][3]);
                slprb[band_index] *= 0.001;   /* vs / 1000 */

                intrb[band_index] = intr[band_index][0]
                    + u*(intr[band_index][2] - intr[band_index][0])
                    + v*(intr[band_index][1] - intr[band_index][0])
                    + u_x_v*(intr[band_index][0] - intr[band_index][1]
                             - intr[band_index][2] + intr[band_index][3]);
                intrb[band_index] *= 0.001;   /* vs / 1000 */
            }

            /* Compute the NDWI variables */
            ndwi_th1 = (andwi[ratio_pix[0]] + 2.0 *
                        sndwi[ratio_pix[0]]) * 0.001;
            ndwi_th2 = (andwi[ratio_pix[0]] - 2.0 *
                        sndwi[ratio_pix[0]]) * 0.001;

            /* Calculate NDWI variables for the band ratios */
            xndwi = ((double) sband[SR_BAND5][nearest_pix] -
                     (double) (sband[SR_BAND7][nearest_pix] * 0.5)) /
                    ((double) sband[SR_BAND5][nearest_pix] +
                     (double) (sband[SR_BAND7][nearest_pix] * 0.5));

            if (xndwi > ndwi_th1)
                xndwi = ndwi_th1;
            if (xndwi < ndwi_th2)
                xndwi = ndwi_th2;

            /* Initialize the band ratios */
            for (band_index = 0; band_index < NSR_BANDS; band_index++)
            {
                erelc[band_index] = -1.0;
                troatm[band_index] = 0.0;
            }

            /* Compute the band ratio */
            erelc[DN_BAND1] = xndwi * slprb[0] + intrb[0];
            erelc[DN_BAND2] = xndwi * slprb[1] + intrb[1];
            erelc[DN_BAND4] = 1.0;
            erelc[DN_BAND7] = xndwi * slprb[2] + intrb[2];

            /* Retrieve the TOA reflectance values for the current pixel */
            troatm[DN_BAND1] = aerob1[nearest_pix];
            troatm[DN_BAND2] = aerob2[nearest_pix];
            troatm[DN_BAND4] = aerob4[nearest_pix];
            troatm[DN_BAND7] = aerob7[nearest_pix];

            /* Retrieve the aerosol information for eps 1.0 */
            iband1 = DN_BAND4;
            iband3 = DN_BAND1;
            eps1 = 1.0;
            iaots = 0;
            subaeroret_new (iband1, iband3, erelc, troatm, tgo_arr,
                xrorayp_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                normext_p0a3_arr, &raot, &residual, &iaots, eps1);

            /* Save the data */
            residual1 = residual;
            sraot1 = raot;

            /* Retrieve the aerosol information for eps 1.75 */
            eps2 = 1.75;
            subaeroret_new (iband1, iband3, erelc, troatm, tgo_arr,
                xrorayp_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                normext_p0a3_arr, &raot, &residual, &iaots, eps2);

            /* Save the data */
            residual2 = residual;

            /* Retrieve the aerosol information for eps 2.5 */
            eps3 = 2.5;
            subaeroret_new (iband1, iband3, erelc, troatm, tgo_arr,
                xrorayp_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                normext_p0a3_arr, &raot, &residual, &iaots, eps3);

            /* Save the data */
            residual3 = residual;
            sraot3 = raot;

            /* Find the eps that minimizes the residual.  This is performed
               by applying a parabolic (quadratic) fit to the three
               (epsilon, residual) pairs found above:
                   r = a\eps^2 + b\eps + c
               The minimum occurs where the first derivative is zero:
                   r' = 2a\eps + b = 0
                   \eps_min = -b/2a

               The a and b coefficients are solved for in the three
               r (residual) equations by eliminating c:
                   r_1 - r_3 = a(\eps_1^2 - \eps_3^2) + b(\eps_1 - \eps_3)
                   r_2 - r_3 = a(\eps_2^2 - \eps_3^2) + b(\eps_2 - \eps_3) */
            xa = (residual1 - residual3)*(eps2 - eps3);
            xb = (residual2 - residual3)*(eps1 - eps3);
            epsmin = 0.5*(xa*(eps2 + eps3) - xb*(eps1 + eps3))/(xa - xb);
            teps[center_pix] = epsmin;

            if (epsmin >= 1.0 && epsmin <= 2.5)
            {
                subaeroret_new (iband1, iband3, erelc, troatm, tgo_arr,
                    xrorayp_arr, roatm_iaMax, roatm_coef, ttatmg_coef,
                    satm_coef, normext_p0a3_arr, &raot, &residual, &iaots,
                    epsmin);
            }
            else if (epsmin < 1.0)
            {
                teps[center_pix] = eps1;
                residual = residual1;
                raot = sraot1;
            }
            else /* epsmin > 2.5 */
            {
                teps[center_pix] = eps3;
                residual = residual3;
                raot = sraot3;
            }

            taero[center_pix] = raot;
            corf = raot / xmus;

            /* Check the model residual.  Corf represents aerosol impact.
               Test the quality of the aerosol inversion. */
            if (residual < (0.015 + 0.005 * corf + 0.10 * troatm[DN_BAND7]))
            {
                /* Test if band 5 makes sense */
                iband = DN_BAND5;
                rotoa = aerob5[nearest_pix];
                atmcorlamb2_new (tgo_arr[iband], xrorayp_arr[iband],
                    aot550nm[roatm_iaMax[iband]], &roatm_coef[iband][0],
                    &ttatmg_coef[iband][0], &satm_coef[iband][0], raot,
                    iband, normext_p0a3_arr[iband], rotoa, &ros5,
                    teps[center_pix]);

                /* Test if band 4 makes sense */
                iband = DN_BAND4;
                rotoa = aerob4[nearest_pix];
                atmcorlamb2_new (tgo_arr[iband], xrorayp_arr[iband],
                    aot550nm[roatm_iaMax[iband]], &roatm_coef[iband][0],
                    &ttatmg_coef[iband][0], &satm_coef[iband][0], raot,
                    iband, normext_p0a3_arr[iband], rotoa, &ros4,
                    teps[center_pix]);

                /* Use the NDVI to validate the reflectance values */
                if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                {
                    /* Clear pixel with valid aerosol retrieval */
                    taero[center_pix] = raot;
                    ipflag[center_pix] |= (1 << IPFLAG_CLEAR);
                }
                else
                {
                    /* Flag as water and use generic values */
                    ipflag[center_pix] |= (1 << IPFLAG_WATER);
                    taero[center_pix] = DEFAULT_AERO;
                    teps[center_pix] = DEFAULT_EPS;
                }
            }
            else
            {
                /* Flag as water and use generic values */
                ipflag[center_pix] |= (1 << IPFLAG_WATER);
                taero[center_pix] = DEFAULT_AERO;
                teps[center_pix] = DEFAULT_EPS;
            }
        }  /* center sample loop */
    }  /* center line loop */

    /* Done with the aerob* arrays */
    free (aerob1);  aerob1 = NULL;
    free (aerob2);  aerob2 = NULL;
    free (aerob4);  aerob4 = NULL;
    free (aerob5);  aerob5 = NULL;
    free (aerob7);  aerob7 = NULL;

    /* Done with the ratiob* arrays */
    free (andwi);  andwi = NULL;
    free (sndwi);  sndwi = NULL;
    free (ratiob1);  ratiob1 = NULL;
    free (ratiob2);  ratiob2 = NULL;
    free (ratiob7);  ratiob7 = NULL;
    free (intratiob1);  intratiob1 = NULL;
    free (intratiob2);  intratiob2 = NULL;
    free (intratiob7);  intratiob7 = NULL;
    free (slpratiob1);  slpratiob1 = NULL;
    free (slpratiob2);  slpratiob2 = NULL;
    free (slpratiob7);  slpratiob7 = NULL;

    /* Done with the DEM, water vapor, and ozone arrays */
    free (dem);  dem = NULL;
    free (wv);  wv = NULL;
    free (oz);  oz = NULL;

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Find the median of the clear aerosols */
    mytime = time(NULL);
    printf ("Computing median of clear pixels in NxN windows %s",
        ctime(&mytime));
    median_aerosol = find_median_aerosol (ipflag, taero, nlines, nsamps);
    if (median_aerosol == 0.0)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("Median aerosol value for clear aerosols is %f\n", median_aerosol);

    /* Fill the cloud, shadow, and water pixels with the median aerosol
       value instead of the default aerosol value */
    mytime = time(NULL);
    printf ("Fill non-clear aerosol values in NxN windows with the median %s",
        ctime(&mytime));
    aerosol_fill_median (ipflag, taero, median_aerosol, nlines, nsamps);

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag2.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols2.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Use the center of the aerosol windows to interpolate the remaining
       pixels in the window */
    mytime = time(NULL);
    printf ("Interpolating the aerosol values in the NxN windows %s",
        ctime(&mytime));
    aerosol_interp (xml_metadata, qaband, ipflag, taero, median_aerosol,
        nlines, nsamps);

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag3.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols3.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Use the center of the aerosol windows to interpolate the teps values
       (angstrom coefficient).  The median value used for filling in clouds and
       water will be the default eps value. */
    mytime = time(NULL);
    printf ("Interpolating the teps values in the NxN windows %s",
        ctime(&mytime));
    aerosol_interp (xml_metadata, qaband, ipflag, teps, DEFAULT_EPS,
        nlines, nsamps);

    /* Perform the second level of atmospheric correction using the aerosols */
    mytime = time(NULL);
    printf ("Performing atmospheric correction.  %s", ctime(&mytime));

    /* 0 .. DN_BAND7 is the same as 0 .. SR_BAND7 here, since the pan band
       isn't spanned */
    for (ib = DN_BAND1; ib <= DN_BAND7; ib++)
    {
        float *band_ptr = sband[ib];   /* convenience pointer */
        int curr_pix;                  /* current pixel being processed */

        printf ("  Band %d\n", ib+1);
#ifdef _OPENMP
        #pragma omp parallel for private (curr_pix) \
        num_threads(get_num_threads())
#endif
        for (curr_pix = 0; curr_pix < nlines*nsamps; curr_pix++)
        {
            float rsurf;        /* surface reflectance */
            float rotoa;        /* top of atmosphere reflectance */
            float roslamb;      /* lambertian surface reflectance */
            float tmpf;         /* temporary floating point value */

            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* If this pixel is cloud, then don't process. taero values
               are generic values anyhow, but TOA values will be returned
               for clouds (not shadows). */
            if (is_cloud (qaband[curr_pix]))
                continue;

            /* Correct all pixels */
            rsurf = band_ptr[curr_pix];
            rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf) +
                     broatm[ib]) * btgo[ib];
            atmcorlamb2_new (tgo_arr[ib], xrorayp_arr[ib],
                             aot550nm[roatm_iaMax[ib]], roatm_coef[ib],
                             ttatmg_coef[ib], satm_coef[ib], taero[curr_pix],
                             ib, normext_p0a3_arr[ib], rotoa, &roslamb,
                             teps[curr_pix]);

            /* If this is the coastal aerosol band then set the aerosol
               bits in the QA band */
            if (ib == DN_BAND1)
            {
                /* Set up aerosol QA bits */
                tmpf = fabsf(rsurf - roslamb);
                if (tmpf <= 0.015)
                {  /* Set the first aerosol bit (low aerosols) */
                    ipflag[curr_pix] |= (1 << AERO1_QA);
                }
                else if (tmpf < 0.03)
                {  /* Set the second aerosol bit (average aerosols) */
                    ipflag[curr_pix] |= (1 << AERO2_QA);
                }
                else
                {  /* Set both aerosol bits (high aerosols) */
                    ipflag[curr_pix] |= (1 << AERO1_QA);
                    ipflag[curr_pix] |= (1 << AERO2_QA);
                }
            }  /* end if this is the coastal aerosol band */

            /* Save the unscaled surface reflectance value.  The data will be
               scaled for output and range checked in convert_output(). */
            band_ptr[curr_pix] = roslamb;
        }  /* pixel loop */
    }  /* end for ib */

    /* Free memory for arrays no longer needed */
    free (taero);
    free (teps);
 
    /* Write the data to the output file */
    mytime = time(NULL);
    printf ("Writing surface reflectance corrected data to the output "
            "files.  %s", ctime(&mytime));

    /* Open the output file */
    sr_output = open_output (xml_metadata, input, OUTPUT_SR);
    if (sr_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through the reflectance bands and write the data */
    for (ib = DN_BAND1; ib <= DN_BAND7; ib++)
    {
        printf ("  Band %d: %s\n", ib+1,
            sr_output->metadata.band[ib].file_name);
        convert_output (sband, ib, out_band, nlines, nsamps, false);
        if (put_output_lines (sr_output, out_band, ib, 0, nlines,
            sizeof (uint16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output data for band %d", ib);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&sr_output->metadata.band[ib],
            &xml_metadata->global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Write the ENVI header */
        strcpy (envi_file, sr_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Append the surface reflectance bands (1-7) to the XML file */
    if (append_metadata (7, sr_output->metadata.band, xml_infile) !=
        SUCCESS)
    {
        sprintf (errmsg, "Appending surface reflectance bands to the "
            "XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the aerosol QA band */
    printf ("  Band %d: %s\n", SR_AEROSOL+1,
            sr_output->metadata.band[SR_AEROSOL].file_name);
    if (put_output_lines (sr_output, ipflag, SR_AEROSOL, 0, nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing aerosol QA output data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Free memory for ipflag data */
    free (ipflag);

    /* Create the ENVI header for the aerosol QA band */
    if (create_envi_struct (&sr_output->metadata.band[SR_AEROSOL],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SR_AEROSOL].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Append the aerosol QA band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SR_AEROSOL],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending aerosol QA band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the output surface reflectance products */
    close_output (sr_output, OUTPUT_SR);
    free_output (sr_output, OUTPUT_SR);

    /* Free the spatial mapping pointer */
    free (space);

    /* Free the data arrays */
    free (rolutt);
    free (transt);
    free (sphalbt);
    free (normext);
    free (tsmax);
    free (tsmin);
    free (nbfic);
    free (nbfi);
    free (ttv);
    free (out_band);

    /* Successful completion */
    mytime = time(NULL);
    printf ("Surface reflectance correction complete.  %s\n", ctime(&mytime));
    return (SUCCESS);
}


/******************************************************************************
MODULE:  is_cloud

PURPOSE:  Determines if the pixel is a cloud (cloud or cirrus cloud).  The
Level-1 QA band is used.  A confidence of high for either of the QA types will
result in the pixel being flagged as cloudy.

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           Pixel is not cloud
true            Pixel is cloud

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool is_cloud
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    /* If the confidence level is high for cloud or cirrus, then flag this as
       a cloud */
    if (level1_qa_cloud_confidence (l1_qa_pix) == L1QA_HIGH_CONF ||
        level1_qa_cirrus_confidence (l1_qa_pix) == L1QA_HIGH_CONF)
        return (true);
    else
        return (false);
}


/******************************************************************************
MODULE:  is_cloud_or_shadow

PURPOSE:  Determines if the pixel is a cloud (cloud, cloud shadow, or cirrus
cloud).  The Level-1 QA band is used.  A confidence of high for any of the
three QA types will result in the pixel being flagged as cloudy.

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           Pixel is not cloud or shadow
true            Pixel is cloud or shadow

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool is_cloud_or_shadow
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    /* If the confidence level is high for cloud, cloud shadow, or cirrus, then
       flag this as a cloud */
    if (level1_qa_cloud_confidence (l1_qa_pix) == L1QA_HIGH_CONF ||
        level1_qa_cloud_shadow_confidence (l1_qa_pix) == L1QA_HIGH_CONF ||
        level1_qa_cirrus_confidence (l1_qa_pix) == L1QA_HIGH_CONF)
        return (true);
    else
        return (false);
}


/******************************************************************************
MODULE:  is_shadow

PURPOSE:  Determines if the pixel is a cloud shadow.  The Level-1 QA band is
used.  A confidence of high for this QA type will result in the pixel being
flagged as a cloud shadow.

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           Pixel is not cloud shadow
true            Pixel is cloud shadow

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool is_shadow
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    /* If the confidence level is high for cloud shadow, then flag this as a
       cloud */
    if (level1_qa_cloud_shadow_confidence (l1_qa_pix) == L1QA_HIGH_CONF)
        return (true);
    else
        return (false);
}
