/*!
\file
\brief Parsing of command-line arguments

\date   Started 1/1/17
\author George
\version\verbatim $Id: cmdline.c 21158 2017-06-10 22:36:08Z karypis $ \endverbatim
*/


#include "kt.h"

/*-------------------------------------------------------------------
 * Command-line options
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"kttype",            1,      0,      CMD_KTTYPE},
  {"iftype",            1,      0,      CMD_IFTYPE},

  {"jbsize",            1,      0,      CMD_JBSIZE},
  {"ibsize",            1,      0,      CMD_IBSIZE},

  {"seed",              1,      0,      CMD_SEED},

  {"dbglvl",            1,      0,      CMD_DBGLVL},
  {"stream",            1,      0,      CMD_STREAM},
  {"help",              0,      0,      CMD_HELP},
  {0,                   0,      0,      0}
};

/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static gk_StringMap_t kttype_options[] = {
  {"baseline1",          KTTYPE_BASELINE1},
  {"baseline2",          KTTYPE_BASELINE2},
  {"baseline3",          KTTYPE_BASELINE3},
  {"baseline4",          KTTYPE_BASELINE4},
  {"baseline5",          KTTYPE_BASELINE5},
  {"baseline6",          KTTYPE_BASELINE6},
  {"baseline7",          KTTYPE_BASELINE7},
  {"baseline8",          KTTYPE_BASELINE8},
  {"baseline9",          KTTYPE_BASELINE9},
  {"baseline10",         KTTYPE_BASELINE10},
  {"baseline5b",         KTTYPE_BASELINE5b},
  {"baseline10b",        KTTYPE_BASELINE10b},
  {"baseline10g",        KTTYPE_BASELINE10g},
  {"baseline5c",         KTTYPE_BASELINE5c},
  {"twopass",            KTTYPE_TWOPASS},
  {"snd",                KTTYPE_SND},
  {"and",                KTTYPE_AND},
  {NULL,                 0}
};

static gk_StringMap_t iftype_options[] = {
  {"tsv",         IFTYPE_TSV},
  {"metis",       IFTYPE_METIS},
  {NULL,                 0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: kt [options] infile [outfile]",
" ",
" ",
" Options",
"  -kttype=text",
"     Specifies the type of k-truss algorithm to use.",
"     Possible values are:",
"        twopass        two-pass parallel [default]",
"        snd            synchronous nucleus decomposition (SND)",
"        and            asynchronous nucleus decomposition (AND)",
"        baseline1      the initial version ",
"        baseline2      the initial version ",
"        baseline3      the initial version ",
"        baseline4      the initial version ",
"        baseline5      the initial version ",
"        baseline5b     the initial version ",
"        baseline5c     the initial version ",
"        baseline6      the initial version ",
"        baseline7      the initial version ",
"        baseline8      the initial version ",
"        baseline9      the initial version ",
"        baseline10     the initial version ",
"        baseline10b    the initial version ",
"        baseline10g    the initial version ",
" ",
"  -iftype=text",
"     Specifies the format of the input file. ",
"     Possible values are:",
"        tsv     tsv format (i, j, v)",
"        metis   Metis format [default]",
" ",
"  -jbsize=int",
"     Specifies the size in words of the J block.",
"     Default value is 1000000.",
" ",
"  -ibsize=int",
"     Specifies the size in words of the I block.",
"     Default value is 1000000.",
" ",
"  -seed=int",
"     Specifies the seed for clustering.",
"     Default value is 1.",
" ",
"  -dbglvl=int",
"     Specifies the level of debugging information to be displayed.",
"     Default value is 0.",
" ",
"  -help",
"     Prints this message.",
""
};



/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
params_t *getcmdline_params(int argc, char *argv[])
{
  gk_idx_t i, j, k;
  int type=0;
  int c, option_index;
  params_t *params;

  params = (params_t *)gk_malloc(sizeof(params_t), "cmdline_parse: params");
  memset(params, 0, sizeof(params_t));

  /* print the command line */
  for (i=0; i<argc; i++)
    printf("%s ", argv[i]);
  printf("\n");

  /* initialize the params data structure */
  params->infile  = NULL;
  params->outfile = NULL;

  params->kttype   = KTTYPE_TWOPASS;
  params->iftype   = IFTYPE_METIS;
  params->jbsize   = 1000000;
  params->ibsize   = 1000000;
  params->dbglvl   = 0;
  params->stream   = 0;

  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_KTTYPE:
        if (gk_optarg) {
          if ((params->kttype = gk_GetStringID(kttype_options, gk_optarg)) == -1)
            errexit("Invalid kttype of %s.\n", gk_optarg);
        }
        break;

      case CMD_IFTYPE:
        if (gk_optarg) {
          if ((params->iftype = gk_GetStringID(iftype_options, gk_optarg)) == -1)
            errexit("Invalid iftype of %s.\n", gk_optarg);
        }
        break;


      case CMD_JBSIZE:
        if (gk_optarg) {
          if ((params->jbsize = atoi(gk_optarg)) < 1)
            errexit("The -jbsize must be greater than 0.\n");
        }
        break;

      case CMD_IBSIZE:
        if (gk_optarg) {
          if ((params->ibsize = atoi(gk_optarg)) < 1)
            errexit("The -ibsize must be greater than 0.\n");
        }
        break;


      case CMD_SEED:
        if (gk_optarg) {
          if ((params->seed = atoi(gk_optarg)) < 1)
            errexit("The -seed must be greater than 0.\n");
        }
        break;

      case CMD_DBGLVL:
        if (gk_optarg) {
          params->dbglvl = atoi(gk_optarg);
          if (params->dbglvl < 0)
            errexit("The -dbglvl must be non-negative.\n");
        }
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(EXIT_SUCCESS);
        break;

      case CMD_STREAM:
        params->stream = atoi(gk_optarg);
        break;

      default:
        printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  /* Get the operation to be performed */
  if (argc-gk_optind == 0 || argc-gk_optind > 2) {
    printf("Missing required parameters.\n  Use %s -help for a summary of the options.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  params->infile = gk_strdup(argv[gk_optind++]);
  if (argc-gk_optind > 0)
    params->outfile = gk_strdup(argv[gk_optind++]);

  return params;
}
