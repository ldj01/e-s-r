#include "aero_interp.h"
#include "quick_select.h"
#include "read_level1_qa.h"
#include "read_level2_qa.h"

/******************************************************************************
MODULE:  aerosol_interp

PURPOSE:  Interpolates the aerosol values throughout the image using the
aerosols that were calculated for each NxN window. Also cleans up the fill
pixels in the ipflag.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_interp
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will fill in the pixels for the remaining, non-center
                          pixels of the window. */
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
)
{
    int i;                 /* looping variable for the bands */
    int line, samp;        /* looping variable for lines and samples */
    int curr_pix;          /* current pixel in 1D arrays of nlines * nsamps */
    int center_line;       /* line for the center of the aerosol window */
    int center_line1;      /* line+1 for the center of the aerosol window */
    int center_lindex;     /* center line array index */
    int center_lindex1;    /* center line 1 array index */
    int center_samp;       /* sample for the center of the aerosol window */
    int center_samp1;      /* sample+1 for the center of the aerosol window */
    int refl_indx = -99;   /* index of band 1 or first band */
    int aero_pix11;        /* pixel location for aerosol window values
                              [lcmg][scmg] */
    int aero_pix12;        /* pixel location for aerosol window values
                              [lcmg][scmg2] */
    int aero_pix21;        /* pixel location for aerosol window values
                              [lcmg2][scmg] */
    int aero_pix22;        /* pixel location for aerosol window values
                              [lcmg2][scmg2] */
    float xaero, yaero;    /* x/y location for aerosol pixel within the overall
                              larger aerosol window grid */
    float aero11;          /* aerosol value at window line, samp */
    float aero12;          /* aerosol value at window line, samp+1 */
    float aero21;          /* aerosol value at window line+1, samp */
    float aero22;          /* aerosol value at window line+1, samp+1 */
    float u, v;            /* line, sample fractional distance from current
                              pixel (weight applied to furthest line, sample) */
    float u_x_v;           /* u * v */
    int aero_window_index_step = AERO_WINDOW*nsamps; /* aerosol window array
                                                        step size */
    float aero_step = 1.0/AERO_WINDOW; /* fraction of window size representing
                                          1 pixel */

    /* Use band 1 band-related metadata for the reflectance information for
       Landsat (Level 1 products).  If band 1 isn't available then just use the
       first band in the XML file. */
    for (i = 0; i < xml_metadata->nbands; i++)
    {
        if (!strcmp (xml_metadata->band[i].name, "b1") &&
            !strncmp (xml_metadata->band[i].product, "L1", 2))
        {
            /* this is the index we'll use for Landsat band info */
            refl_indx = i;
        }
    }
    if (refl_indx == -99)
        refl_indx = 0;

    /* Interpolate the aerosol data for each pixel location */
    for (line = 0, curr_pix = 0; line < nlines; line++)
    {
        /* Determine the line of the representative center pixel in the
           aerosol NxN window array */
        center_line = (int)(line*aero_step)*AERO_WINDOW + HALF_AERO_WINDOW;
        center_lindex = center_line*nsamps;

        /* Determine fractional location of this line in the aerosol window.
           Negative values are at the top of the window. */
        yaero = (line - center_line)*aero_step;
        u = yaero - (int) yaero;

        /* Determine if this pixel is closest to the line below or the line
           above. If the fractional value is in the top part of the aerosol
           window, then use the line above.  Otherwise use the line below. */
        if (u < 0.0)
        {
            center_line1 = center_line - AERO_WINDOW;
            center_lindex1 = center_lindex - aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 < 0)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }
        else
        {
            center_line1 = center_line + AERO_WINDOW;
            center_lindex1 = center_lindex + aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 >= nlines - 1)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }

        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* If this pixel is cloud or shadow, then don't process. Use
               median aerosol values.  Flag them separately. */
            else if (is_cloud (qaband[curr_pix]))
            {
                taero[curr_pix] = median_aero;
                ipflag[curr_pix] = (1 << IPFLAG_CLOUD);
                continue;
            }
            else if (is_shadow (qaband[curr_pix]))
            {
                taero[curr_pix] = median_aero;
                ipflag[curr_pix] = (1 << IPFLAG_SHADOW);
                continue;
            }

            /* If this pixel is water, then don't process. Use default aerosol
               values. */
            else if (level1_qa_is_water (qaband[curr_pix]))
            {
                taero[curr_pix] = median_aero;
                ipflag[curr_pix] = (1 << IPFLAG_WATER);
                continue;
            }

            /* Determine the sample of the representative center pixel in the
               aerosol NxN window array */
            center_samp = (int)(samp*aero_step)*AERO_WINDOW + HALF_AERO_WINDOW;

            /* If the current line, sample are the center line, sample, then
               skip to the next pixel.  We already have the aerosol value. */
            if (samp == center_samp && line == center_line)
                continue;

            /* Determine fractional location of this sample in the aerosol
               window.  Negative values are at the left of the window. */
            xaero = (samp - center_samp)*aero_step;
            v = xaero - (int) xaero;

            /* Determine if this pixel is closest to the sample to the left or
               the sample to the right.  If the fractional value is on the left
               side of the aerosol window, then use the sample to the left.
               Otherwise use the sample to the right. */
            if (v < 0.0)
            {
                center_samp1 = center_samp - AERO_WINDOW;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 < 0)
                    center_samp1 = center_samp;
            }
            else
            {
                center_samp1 = center_samp + AERO_WINDOW;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 >= nsamps-1)
                    center_samp1 = center_samp;
            }

            /* Determine the four aerosol window pixels to be used for
               interpolating the current pixel */
            aero_pix11 = center_lindex + center_samp;
            aero_pix12 = center_lindex + center_samp1;
            aero_pix21 = center_lindex1 + center_samp;
            aero_pix22 = center_lindex1 + center_samp1;

            /* Get the aerosol values */
            aero11 = taero[aero_pix11];
            aero12 = taero[aero_pix12];
            aero21 = taero[aero_pix21];
            aero22 = taero[aero_pix22];

            /* From here make the fractional distance positive, regardless of
               where it is in the window. */
            u = fabsf(u);
            v = fabsf(v);
            u_x_v = u * v;

            /* Interpolate the aerosol */
            taero[curr_pix] = aero11
                            + u*(aero21 - aero11)
                            + v*(aero12 - aero11)
                            + u_x_v*(aero11 - aero12 - aero21 + aero22);

            /* Set the aerosol to window interpolated. Clear anything else. */
            ipflag[curr_pix] = (1 << IPFLAG_INTERP_WINDOW);

            /* If any of the window pixels used in the interpolation were
               water pixels, then mask this pixel with water (in addition to
               the interpolation bit already set) */
            if (lasrc_qa_is_water (ipflag[aero_pix11]) ||
                lasrc_qa_is_water (ipflag[aero_pix12]) ||
                lasrc_qa_is_water (ipflag[aero_pix21]) ||
                lasrc_qa_is_water (ipflag[aero_pix22]))
                ipflag[curr_pix] |= (1 << IPFLAG_WATER);
        }  /* end for samp */
    }  /* end for line */

    /* Clean up the ipflag in the center of the NxN windows, for the fill
       pixels. If an NxN window is a mixture of fill and non-fill, the center
       of the window can be flagged as fill and some other QA based on the
       other pixels in that window. At the end, we want fill to be fill. */
    for (curr_pix = 0; curr_pix < nlines * nsamps; curr_pix++)
    {
        if (level1_qa_is_fill (qaband[curr_pix]))
            ipflag[curr_pix] = (1 << IPFLAG_FILL);
    }
}


/******************************************************************************
MODULE:  aerosol_fill_median

PURPOSE:  Changes the aerosol window values for water, cloud, and shadow to
be the median aerosol value of the clear pixels.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_fill_median
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will interpolate/average the pixels of the windows
                          that failed the aerosol inversion (using ipflag) */
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
)
{
    int line, samp;       /* looping variable for lines and samples */
    int curr_pix;         /* current pixel in 1D arrays of nlines * nsamps */

    /* Loop through the center of the NxN window pixels */
    for (line = HALF_AERO_WINDOW; line < nlines; line += AERO_WINDOW)
    {
        curr_pix = line * nsamps + HALF_AERO_WINDOW;
        for (samp = HALF_AERO_WINDOW; samp < nsamps;
             samp += AERO_WINDOW, curr_pix += AERO_WINDOW)
        {
            /* Find cloud, shadow, and water pixels and reset the default
               aerosol value to that of the median aerosol value */
            if (lasrc_qa_is_cloud_cirrus (ipflag[curr_pix]) ||
                lasrc_qa_is_cloud_shadow (ipflag[curr_pix]) ||
                lasrc_qa_is_water (ipflag[curr_pix]))
            {
                taero[curr_pix] = median_aero;
            }
        }
    }
}


/******************************************************************************
MODULE:  find_median_aerosol

PURPOSE:  Finds the median aerosol value for the valid land aerosols.

RETURN VALUE:
Type = float
Value           Description
-----           -----------
zero            Error allocating memory for the aerosol array
non-zero        Median aerosol value of the array

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
float find_median_aerosol
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the center of the aerosol
                             windows. */
    float *taero,      /* I: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the center of the aerosol windows */
    int nlines,        /* I: number of lines in taero band */
    int nsamps         /* I: number of samps in taero band */
)
{
    char errmsg[STR_SIZE];                         /* error message */
    char FUNC_NAME[] = "find_median_aerosol";      /* function name */
    int line, samp;       /* looping variable for lines and samples */
    int curr_pix;         /* current pixel in 1D arrays of nlines * nsamps */
    int nbclrpix;         /* number of clear aerosol pixels in this array */
    int nwindows;         /* number of NxN windows in the image */
    float median;         /* median clear aerosol value */
    float *aero = NULL;   /* array of the clear aerosol values */

    /* Determine how many NxN windows there are in this array of data */
    nwindows = ceil ((float) nlines / AERO_WINDOW) *
               ceil ((float) nsamps / AERO_WINDOW);

    /* Allocate memory for the aerosols in each window */
    aero = calloc (nwindows, sizeof (float));
    if (aero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for clear aerosol array");
        error_handler (true, FUNC_NAME, errmsg);
        return (0.0);
    }

    /* Loop through the NxN center window values and write the clear aerosol
       values to the aerosol array for determining the median */
    nbclrpix = 0;
    for (line = HALF_AERO_WINDOW; line < nlines; line += AERO_WINDOW)
    {
        curr_pix = line * nsamps + HALF_AERO_WINDOW;
        for (samp = HALF_AERO_WINDOW; samp < nsamps;
             samp += AERO_WINDOW, curr_pix += AERO_WINDOW)
        {
            /* Process clear aerosols */
            if (lasrc_qa_is_valid_aerosol_retrieval (ipflag[curr_pix]))
            {
                aero[nbclrpix] = taero[curr_pix];
                nbclrpix++;
            }  /* if pixel is clear */
        }  /* for samp */
    }  /* for line */

    /* If no clear aerosols were available, then just return a default value */
    if (nbclrpix == 0)
        median = DEFAULT_AERO;
    else
    {
        /* Get the median of the clear pixels */
        median = quick_select (aero, nbclrpix);
    }

    /* Free memory */
    free (aero);

    /* Successful completion */
    return (median);
}
