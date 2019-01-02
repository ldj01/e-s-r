#include <sys/stat.h>
#include <unistd.h>
#include "lasrc.h"

/******************************************************************************
MODULE:  lasrc (Landsat Surface Reflectance Code - LaSRC)

PURPOSE:  Computes the surface reflectance values for the Landsat 8 products.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the surface reflectance
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Bands 1-7 are corrected to surface reflectance.  Band 8 (pand band) is not
   processed.  Band 9 (cirrus band) is corrected to TOA reflectance.  Bands
   10 and 11 are corrected to brightness temperature.
2. SDstart and SDreaddata have minor memory leaks.  Ultimately both call
   HAregister_atom which makes a malloc call and the memory is never freed.
3. Conversion algorithms for TOA reflectance and TOA brightness temperature are
   available from
   http://landsat.usgs.gov/Landsat8_Using_Product.php
4. The TOA and SR corrections utilize solar and view angles.  Previous versions
   simply corrected every pixel for the solar angles at the scene center.
   However as of version 1.0.0, the solar and view angles are read from the
   per-pixel angle bands and better reflect the angle at that location within
   the scene.  These solar/view angles are scaled int16s and therefore need to
   be unscaled (multiply by 0.01) before using.  They are in units of degrees.
5. The solar angles from band to band are fairly stable, thus a single array of
   per-pixel angles will be used.  The array can be from a "representative
   band" for the reflectance bands or it can be an average of the reflectance
   bands. We will use band 4 as the representative band for these per-pixel
   angle values.
6. The view/observation angles from band to band are a bit more unstable
   (particularly the view azimuth) near nadir.  A single array of per-pixel
   angles will still be used, however the information near nadir may be slightly
   affected.  We will use band 4 as the representative band for these per-pixel
   angle values.
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;            /* verbose flag for printing messages */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char envi_file[STR_SIZE];/* ENVI filename */
    char *aux_path = NULL;   /* path for Landsat auxiliary data */
    char *xml_infile = NULL; /* input XML filename */
    char *aux_infile = NULL; /* input auxiliary filename for water vapor
                                and ozone*/
    char *cptr = NULL;       /* pointer to the file extension */
    char aux_year[5];        /* string to contain the year of auxiliary file */

    int retval;              /* return status */
    int ib;                  /* looping variable for input bands */
    int i;                   /* looping variables */
    Input_t *input = NULL;       /* input structure for the Landsat product */
    Output_t *toa_output = NULL; /* output structure and metadata for the TOA
                                    product */
    Output_t *radsat_output = NULL; /* output structure and metadata for the
                                       RADSAT product */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global meta */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    struct stat statbuf;      /* buffer for the file stat function */

    int16 *sza = NULL;  /* per-pixel solar zenith angles, nlines x nsamps */
    int16 *saa = NULL;  /* per-pixel solar azimuth angles, nlines x nsamps */
    int16 *vza = NULL;  /* per-pixel view zenith angles, nlines x nsamps */
    int16 *vaa = NULL;  /* per-pixel view azimuth angles, nlines x nsamps */
    float **sband = NULL;     /* output surface reflectance and brightness
                                 temp bands, qa band is separate as a uint16 */
    uint16 *out_band = NULL;  /* scaled output */
    uint16 *qaband = NULL;    /* QA band for the input image, nlines x nsamps */
    uint16 *radsat = NULL;    /* QA band for radiometric saturation of the
                                 Level-1 product, nlines x nsamps */

    float xts;           /* scene center solar zenith angle (deg) */
    float xmus;          /* cosine of solar zenith angle */
    bool process_sr = true;  /* this is set to false if the solar zenith
                                is too large and the surface reflectance
                                cannot be calculated or if the user specifies
                                that surface reflectance processing will not
                                be completed and only TOA processing will be
                                done */
    bool write_toa = false;  /* this is set to true if the user specifies
                                TOA products should be output for delivery */
    float pixsize;      /* pixel size for the reflectance bands */
    int nlines, nsamps; /* number of lines and samples in the reflectance and
                           thermal bands */

    /* The following arguments are all names of the LUTs. The first five are
       all tables of coefficients generated by the 6S software and provided
       as input to this application. */
    char anglehdf[STR_SIZE];  /* angle HDF filename */
    char intrefnm[STR_SIZE];  /* intrinsic reflectance filename */
    char transmnm[STR_SIZE];  /* transmission filename */
    char spheranm[STR_SIZE];  /* spherical albedo filename */
    char cmgdemnm[STR_SIZE];  /* climate modeling grid DEM filename */
    char rationm[STR_SIZE];   /* ratio averages filename ("ratio map" used by
                                 the aerosol retrieval algorithm) */
    char auxnm[STR_SIZE];     /* auxiliary filename for ozone and water vapor*/

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &aux_infile, &process_sr,
        &write_toa, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    printf ("Starting TOA and surface reflectance processing ...\n");

    /* Provide user information if verbose is turned on */
    if (verbose)
    {
        printf ("  XML input file: %s\n", xml_infile);
        printf ("  AUX input file: %s\n", aux_infile);
        if (!process_sr)
        {
            printf ("    **Surface reflectance corrections will not be "
                "completed.  Only top of atmosphere corrections will be "
                "completed.\n");
        }
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Open the reflectance product, set up the input data structure, and
       allocate memory for the data buffers */
    input = open_input (&xml_metadata, process_sr);
    if (input == (Input_t *) NULL)
    {
        sprintf (errmsg, "Error opening/reading the input DN data: %s",
            xml_infile);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    gmeta = &xml_metadata.global;

    /* Output some information from the input files if verbose */
    if (verbose)
    {
        printf ("  Nband: %d\n", input->nband);
        printf ("  Number of lines/samples: %d/%d\n", input->size.nlines,
            input->size.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size.pixsize[0],
            input->size.pixsize[1]);

        printf ("  Nband thermal: %d\n", input->nband_th);
        printf ("  Number of thermal lines/samples: %d/%d\n",
            input->size_th.nlines, input->size_th.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size_th.pixsize[0],
            input->size_th.pixsize[1]);

        printf ("  Nband QA: %d\n", input->nband_qa);
        printf ("  Number of qa lines/samples: %d/%d\n",
            input->size_qa.nlines, input->size_qa.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size_qa.pixsize[0],
            input->size_qa.pixsize[1]);

        printf ("  Fill value: %d\n", input->meta.fill);
        printf ("  Solar zenith: %f\n", xml_metadata.global.solar_zenith);
        printf ("  Solar azimuth: %f\n", xml_metadata.global.solar_azimuth);
    }

    /* Pull the needed metadata from the XML file and input structure */
    xts = gmeta->solar_zenith;
    pixsize = (float) input->size.pixsize[0];
    xmus = cos (xts * DEG2RAD);
    nlines = input->size.nlines;
    nsamps = input->size.nsamps;

    /* If this is OLI-only data, then surface reflectance can not be
       processed */
    if (input->meta.inst == INST_OLI && process_sr)
    {
        sprintf (errmsg, "This is an OLI-only scene vs. an OLI-TIRS scene. "
            "Corrections must be limited to top-of-atmosphere reflectance and "
            "brightness temperature corrections. Use the --process_sr=false "
            "command-line argument to process. (oli-only cannot be corrected "
            "to surface reflectance)");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* The surface reflectance algorithm cannot be implemented for solar
       zenith angles greater than 76 degrees.  Need to flag if the current
       scene falls into that category. */
    if (xts > 76.0 && process_sr)
    {
        sprintf (errmsg, "Solar zenith angle is too large to allow for surface "
            "reflectance processing.  Corrections must be limited to top-of-"
            "atmosphere reflectance and brightness temperature corrections. "
            "Use the --process_sr=false command-line argument. "
            "(solar zenith angle out of range)");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for all the data arrays */
    if (verbose)
        printf ("Allocating memory for the data arrays ...\n");
    retval = memory_allocation_main (nlines, nsamps, &sza, &saa, &vza, &vaa,
        &qaband, &radsat, &sband, &out_band);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        sprintf (errmsg, "Error allocating memory for the data arrays from "
            "the main application.");
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the QA band */
    if (get_input_qa_lines (input, 0, 0, nlines, qaband) != SUCCESS)
    {
        sprintf (errmsg, "Reading QA band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the scaled solar and view azimuth/zenith per pixel angle bands
       which are in degrees */
    if (get_input_ppa_lines (input, 0, nlines, sza, saa, vza, vaa) != SUCCESS)
    {
        sprintf (errmsg, "Reading per-pixel solar and view angle bands");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Get the L8 auxiliary directory and the full pathname of the auxiliary
       files to be read if processing surface reflectance */
    if (process_sr)
    {
        /* Get the path for the auxiliary products from the L8_AUX_DIR
           environment variable.  If it isn't defined, then assume the products
           are in the local directory. */
        aux_path = getenv ("L8_AUX_DIR");
        if (aux_path == NULL)
        {
            aux_path = ".";
            sprintf (errmsg, "L8_AUX_DIR environment variable isn't defined. "
                "It is assumed the auxiliary products will be available from "
                "the local directory.");
            error_handler (false, FUNC_NAME, errmsg);
        }

        /* Grab the year of the auxiliary input file to be used for the correct
           location of the auxiliary file in the auxiliary directory */
        strncpy (aux_year, &aux_infile[5], 4);
        aux_year[4] = '\0';

        /* Set up the look-up table files and make sure they exist */
        sprintf (anglehdf, "%s/LDCMLUT/ANGLE_NEW.hdf", aux_path);
        sprintf (intrefnm, "%s/LDCMLUT/RES_LUT_V3.0-URBANCLEAN-V2.0.hdf",
            aux_path);
        sprintf (transmnm, "%s/LDCMLUT/TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
            aux_path);
        sprintf (spheranm, "%s/LDCMLUT/AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
            aux_path);
        sprintf (cmgdemnm, "%s/CMGDEM.hdf", aux_path);
        sprintf (rationm, "%s/ratiomapndwiexp.hdf", aux_path);
        sprintf (auxnm, "%s/LADS/%s/%s", aux_path, aux_year, aux_infile);

        if (stat (anglehdf, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find anglehdf data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", anglehdf);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (intrefnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find intrefnm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", intrefnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (transmnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find transmnm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", transmnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (spheranm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find spheranm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", spheranm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (cmgdemnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find cmgdemnm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", cmgdemnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (rationm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find rationm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", rationm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (auxnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find auxnm data file: %s\n  Check "
                "L8_AUX_DIR environment variable.", auxnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Compute the TOA reflectance and TOA brightness temp */
    printf ("Calculating TOA reflectance and TOA brightness temps...\n");
    retval = compute_toa_refl (input, &xml_metadata, qaband, nlines, nsamps,
        gmeta->instrument, sza, sband, radsat);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error computing TOA reflectance and TOA brightness "
            "temperatures.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the TOA output file, and set up the bands according to whether
       the TOA reflectance bands will be written. */
    toa_output = open_output (&xml_metadata, input, OUTPUT_TOA);
    if (toa_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("Writing TOA reflectance corrected data to the output files ...\n");

    /* If we are writing the TOA data, do so now for bands 1-7.  This will
       occur if the user specified TOA to be written or if the surface
       reflectance processing will not be completed. */
    if (write_toa || !process_sr)
    {
        for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
        {
            printf ("  Band %d: %s\n", ib+1,
                toa_output->metadata.band[ib].file_name);
            convert_output (sband, ib, out_band, nlines, nsamps, false);
            if (put_output_lines (toa_output, out_band, ib, 0, nlines,
                sizeof (uint16)) != SUCCESS)
            {
                sprintf (errmsg, "Writing output TOA data for band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Create the ENVI header file this band */
            if (create_envi_struct (&toa_output->metadata.band[ib],
                &xml_metadata.global, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Creating ENVI header structure.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
      
            /* Write the ENVI header */
            strcpy (envi_file, toa_output->metadata.band[ib].file_name);
            cptr = strchr (envi_file, '.');
            strcpy (cptr, ".hdr");
            if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Writing ENVI header file.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Append the TOA reflectance bands, bands 1-7, to the XML file */
        if (append_metadata (7, toa_output->metadata.band, xml_infile) !=
            SUCCESS)
        {
            sprintf (errmsg, "Appending TOA reflectance bands to XML file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Write bands 9-11 (cirrus and thermals), which don't get any further
       processing. */
    for (ib = SR_BAND9; ib <= SR_BAND11; ib++)
    {
        /* If processing OLI-only, then bands 10 and 11 don't exist */
        if (!strcmp (gmeta->instrument, "OLI") &&
            (ib == SR_BAND10 || ib == SR_BAND11))
            continue;
        
        printf ("  Band %d: %s\n", ib+2,
            toa_output->metadata.band[ib].file_name);

        if (ib == SR_BAND9)
            convert_output (sband, ib, out_band, nlines, nsamps, false);
        else 
            convert_output (sband, ib, out_band, nlines, nsamps, true); 

        if (put_output_lines (toa_output, out_band, ib, 0, nlines,
            sizeof (uint16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output TOA data for band %d", ib+2);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&toa_output->metadata.band[ib],
            &xml_metadata.global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
  
        /* Write the ENVI header */
        strcpy (envi_file, toa_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Append the TOA cirrus/thermal band to the XML file */
        if (append_metadata (1, &toa_output->metadata.band[ib], xml_infile) !=
            SUCCESS)
        {
            sprintf (errmsg, "Appending TOA cirrus/thermal band to XML file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the output TOA products, cleanup bands, and free the memory */
    close_output (toa_output, OUTPUT_TOA);
    if (process_sr && !write_toa)
    {
        /* Remove the TOA bands 1-7 that were created by the open routine,
           since they aren't actually used */
        for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
            unlink (toa_output->metadata.band[ib].file_name);
    }
    free_output (toa_output, OUTPUT_TOA);

    /* Open the RADSAT output file */
    radsat_output = open_output (&xml_metadata, input, OUTPUT_RADSAT);
    if (radsat_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("Writing RADSAT data to the output files ...\n");

    /* Write the RADSAT band */
    if (put_output_lines (radsat_output, radsat, 0, 0, nlines,
        sizeof (uint16)) != SUCCESS)
    {
        sprintf (errmsg, "Writing output RADSAT data");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Create the ENVI header file this band */
    if (create_envi_struct (&radsat_output->metadata.band[SR_RADSAT],
        &xml_metadata.global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
      
    /* Write the ENVI header */
    strcpy (envi_file, radsat_output->metadata.band[SR_RADSAT].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Append the RADSAT band to the XML file */
    if (append_metadata (1, &radsat_output->metadata.band[SR_RADSAT],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending the RADSAT band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close the output radsat product, cleanup bands, and free the memory */
    close_output (radsat_output, OUTPUT_RADSAT);
    free_output (radsat_output, OUTPUT_RADSAT);
    free (radsat);

    /* Only continue with the surface reflectance corrections if SR processing
       has been requested and is possible due to the solar zenith angle */
    if (process_sr)
    {
        /* Perform atmospheric correction for the reflectance bands and write
           the data to the SR output file */
        printf ("Performing atmospheric corrections for each reflectance "
            "band ...\n");
        retval = compute_sr_refl (input, &xml_metadata, xml_infile, qaband,
            nlines, nsamps, pixsize, sband, sza, saa, vza, vaa, xts, xmus,
            anglehdf, intrefnm, transmnm, spheranm, cmgdemnm, rationm, auxnm);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Error computing surface reflectance");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end if process_sr */
  
    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the input product */
    if (verbose)
       printf ("Closing input/output and freeing pointers ...\n");
    close_input (input);
    free_input (input);

    /* Free the filename pointers */
    free (xml_infile);
    free (aux_infile);

    /* Free memory for band data */
    free (sza);
    free (saa);
    free (vza);
    free (vaa);
    free (qaband);
    free (out_band);
    for (i = 0; i < NBAND_TTL_OUT-1; i++)
        free (sband[i]);
    free (sband);

    /* Indicate successful completion of processing */
    printf ("Surface reflectance processing complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void usage ()
{
    printf ("LaSRC (Landsat Surface Reflectance Code) computes the surface "
            "reflectance values for the input Landsat 8 DN products.  Surface "
            "reflectance correction and/or top of atmosphere correction is "
            "applied and written for bands 1-7.  Top of atmosphere and "
            "corrections are applied and written for bands 9 (cirrus), "
            "10 (thermal), and 11 (thermal).  Surface reflectance corrections "
            "are available for OLI_TIRS products.  OLI-only scenes are "
            "corrected up through TOA and not surface reflectance.\n\n");
    printf ("usage: lasrc "
            "--xml=input_xml_filename "
            "--aux=input_auxiliary_filename "
            "--process_sr=true:false --write_toa [--verbose] [--version] "
            "[--scale_refl=<X.X>] [--offset_refl=<X.X>] "
            "[--scale_therm=<X.X>] [--offset_therm=<X.X>] \n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");
    printf ("    -aux: name of the input auxiliary file containing ozone "
            "and water vapor for the scene date.  The file is expected to "
            "live in the $L8_AUX_DIR/LADS directory or in the local "
            "directory.\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -process_sr: the default is to process surface reflectance, "
            "however if this flag is set to false then only the TOA "
            "reflectance processing and brightness temperature will be "
            "done.\n");
    printf ("    -write_toa: the intermediate TOA reflectance products "
            "for bands 1-7 are written to the output file\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("    -version: print the LaSRC version. When this parameter is "
            "used, none of the other parameters are used or required.\n");
    printf ("    -scale_refl: Scaling value for reflective bands.\n");
    printf ("    -offset_refl: Offset value for reflective bands.\n");
    printf ("    -scale_therm: Scaling value for thermal bands.\n");
    printf ("    -offset_therm: Offset value for thermal bands.\n");

    printf ("\nlasrc --help will print the usage statement\n");
    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_01_T1.xml "
            "--aux=L8ANC2013181.hdf_fused --verbose\n");
    printf ("   ==> Writes bands 9-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n\n");

    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_01_T1.xml "
            "--aux=L8ANC2013181.hdf_fused --write_toa --verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n");

    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_01_T1.xml "
            "--aux=L8ANC2013181.hdf_fused --process_sr=false --verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Surface reflectance corrections are not applied.\n");
}


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
bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is rightmost bit) */
)
{
    /* Take 2 ** n, then AND that result with the byte value */
    return (byte_val & (1 << n));
}

