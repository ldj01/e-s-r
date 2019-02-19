#include <getopt.h>
#include <stdlib.h>
#include "lasrc.h"

static double scale_refl;    /* scale for reflective bands */
static double offset_refl;   /* add offset for reflective bands */
static double scale_therm;   /* scale for thermal bands */
static double offset_therm;  /* add offset for thermal bands */
static double mult_refl;     /* output reflective scale factor */
static double mult_therm;    /* output thermal scale factor */
static int num_threads;      /* number of threads */

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. The input files should be character a pointer set to NULL on input. Memory
     for these pointers is allocated by this routine. The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
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
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    static int verbose_flag=0;       /* verbose flag */
    static int write_toa_flag=0;     /* write TOA flag */
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "get_args";   /* function name */
    static int version_flag=0;       /* flag to print version number instead
                                        of processing */
    static struct option long_options[] =
    {
        {"verbose", no_argument, &verbose_flag, 1},
        {"write_toa", no_argument, &write_toa_flag, 1},
        {"xml", required_argument, 0, 'i'},
        {"aux", required_argument, 0, 'a'},
        {"process_sr", required_argument, 0, 'p'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, &version_flag, 1},
        {"offset_refl", required_argument, 0, 'm'},
        {"offset_therm", required_argument, 0, 'n'},
        {"scale_refl", required_argument, 0, 'r'},
        {"scale_therm", required_argument, 0, 't'},
        {"num_threads", required_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    /* Initialize the flags to false */
    *verbose = false;
    *write_toa = false;
    *process_sr = true;    /* default is to process SR products */

    /* Initialize to default value */
    scale_refl = SCALE_FACTOR;
    scale_therm = SCALE_FACTOR_TH;
    mult_refl = 1 / SCALE_FACTOR;
    mult_therm = 1 / SCALE_FACTOR_TH;
    offset_refl = OFFSET_REFL;
    offset_therm = OFFSET_THERM;
    num_threads = 1;

    /* Loop through all the cmd-line options */
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 'i':  /* XML input file */
                *xml_infile = strdup (optarg);
                break;
     
            case 'a':  /* auxiliary input file */
                *aux_infile = strdup (optarg);
                break;
     
            case 'p':  /* process SR products */
                if (!strcmp (optarg, "true"))
                    *process_sr = true;
                else if (!strcmp (optarg, "false"))
                    *process_sr = false;
                else
                {
                    sprintf (errmsg, "Unknown value for process_sr: %s",
                        optarg);
                    error_handler (true, FUNC_NAME, errmsg);
                    usage ();
                    return (ERROR);
                }
                break;
     
            case 'm':
                offset_refl = atof(optarg);
                break;

            case 'n':
                offset_therm = atof(optarg);
                break;

            case 'r':
                scale_refl = atof(optarg);
                mult_refl = 1 / scale_refl;
                break;

            case 't':
                scale_therm = atof(optarg);
                mult_therm = 1 / scale_therm;
                break;

            case 'v': /* number of threads */
                num_threads = atoi (optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Print version number instead of processing */
    if (version_flag)
    {
        printf("%s\n", SR_VERSION);
        exit(EXIT_SUCCESS);
    }

    /* Make sure the XML file was specified */
    if (*xml_infile == NULL)
    {
        sprintf (errmsg, "Input XML file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Make sure the auxiliary file was specified */
    if (*aux_infile == NULL)
    {
        sprintf (errmsg, "Input auxiliary file for water vapor and ozone is "
            "a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Check the flags */
    if (verbose_flag)
        *verbose = true;
    if (write_toa_flag)
        *write_toa = true;

    return (SUCCESS);
}

/* Value retrieval functions. */
double get_scale_refl(void)
{
    return scale_refl;
}

double get_scale_therm(void)
{
    return scale_therm;
}

double get_offset_refl(void)
{
    return offset_refl;
}

double get_offset_therm(void)
{
    return offset_therm;
}

double get_mult_refl(void)
{
    return mult_refl;
}

double get_mult_therm(void)
{
    return mult_therm;
}

int get_num_threads(void)
{
    return num_threads;
}
