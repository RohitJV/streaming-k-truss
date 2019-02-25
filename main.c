/*!
\file
\brief The entry point of the triangle counting code
\date Started 5/10/2017
\author George
\version\verbatim $Id: cmdline.c 20946 2017-05-10 23:12:48Z karypis $ \endverbatim
*/

#include "kt.h"

void cleanMem(vault_t *vault, params_t *params) {
  gk_free((void **)&params->infile, LTERM);
  gk_free((void **)&params->outfile, LTERM);
  gk_free((void **)&params, LTERM);

  gk_free((void **)&vault->graph->xadj, LTERM);
  gk_free((void **)&vault->graph->adjncy, LTERM);
  gk_free((void **)&vault->graph, LTERM);

  gk_free((void **)&vault->ugraph->xadj, LTERM);
  gk_free((void **)&vault->ugraph->adjncy, LTERM);
  gk_free((void **)&vault->ugraph->iadjwgt, LTERM);
  gk_free((void **)&vault->ugraph, LTERM);

  gk_free((void **)&vault->lgraph->xadj, LTERM);
  gk_free((void **)&vault->lgraph->adjncy, LTERM);
  gk_free((void **)&vault->lgraph, LTERM);

  gk_free((void **)&vault->perm, LTERM);
  gk_free((void **)&vault->iperm, LTERM);
  gk_free((void **)&vault->ktedges, LTERM);

  gk_free((void **)&vault, LTERM);
}

/*************************************************************************
* The entry point
**************************************************************************/
int main(int argc, char *argv[])
{
  int64_t ntriangles=0;
  params_t *params;
  vault_t *vault;
  params = getcmdline_params(argc, argv);
  double totalRuntime = 0.0;
  setvbuf(stdout, NULL, _IOLBF, 0);

  vault = loadData(params);

  printf("\n-----------------\n");
  printf("infile: %s\n", params->infile);
  printf("outfile: %s\n", (params->outfile ? params->outfile : "NULL"));
  printf("#nvtxs: %d\n", vault->graph->nvtxs);
  printf("kttype: %s\n", kttypenames[params->kttype]);
  printf("jbsize: %d, ibsize: %d\n", params->jbsize, params->ibsize);
  printf("\n");

  /* perform various initializations */
  srand(params->seed);

  gk_startwctimer(vault->timer_global);

  /* perform the k-truss decomposition */
  switch (params->kttype) {
    case KTTYPE_BASELINE1:
      ntriangles = kt_Baseline1(params, vault);
      break;

    case KTTYPE_BASELINE2:
      ntriangles = kt_Baseline2(params, vault);
      break;

    case KTTYPE_BASELINE3:
      ntriangles = kt_Baseline3(params, vault);
      break;

    case KTTYPE_BASELINE4:
      ntriangles = kt_Baseline4(params, vault);
      break;

    case KTTYPE_BASELINE5:
      ntriangles = kt_Baseline5(params, vault);
      break;

    case KTTYPE_BASELINE5b:
      ntriangles = kt_Baseline5b(params, vault);
      break;

    case KTTYPE_BASELINE5c:
      ntriangles = kt_Baseline5c(params, vault);
      break;

    case KTTYPE_BASELINE6:
      ntriangles = kt_Baseline6(params, vault);
      break;

    case KTTYPE_BASELINE7:
      ntriangles = kt_Baseline7(params, vault);
      break;

    case KTTYPE_BASELINE8:
      ntriangles = kt_Baseline8(params, vault);
      break;

#ifdef XXX
    case KTTYPE_BASELINE9:
      ntriangles = kt_Baseline9(params, vault);
      break;
#endif

    case KTTYPE_BASELINE10:
      ntriangles = kt_Baseline10(params, vault);
      break;

    case KTTYPE_BASELINE10b:
      ntriangles = kt_Baseline10b(params, vault);
      break;

    case KTTYPE_BASELINE10g:
      ntriangles = kt_Baseline10g(params, vault);
      break;

    case KTTYPE_TWOPASS:
      ntriangles = kt_twopass5b(params, vault);
      break;

    case KTTYPE_SND:
      ntriangles = kt_snd(params, vault);
      break;

    case KTTYPE_AND:
      ntriangles = kt_and(params, vault);
      break;

    default:
      errexit("Unknown kttype of %d\n", params->kttype);
  }

  gk_stopwctimer(vault->timer_global);

  if (vault->ktedges && params->dbglvl&8)
    kt_CheckKTrussDecomposition(params, vault);

  writeResults(params, vault);

  printf("#triangles: %"PRId64"; rate: %.4lf MT/sec\n", ntriangles,
      ((double)ntriangles)/((double)1e6*gk_getwctimer(vault->timer_2)));
  printf("ktmax (k-2): %d\n", vault->ktmax);
  printf("----------------------------\n");
  if (gk_getwctimer(vault->timer_1) > 0)
    printf("   timer_1: %7.2lfs\n", gk_getwctimer(vault->timer_1));
  if (gk_getwctimer(vault->timer_2) > 0)
    printf("   timer_2: %7.2lfs\n", gk_getwctimer(vault->timer_2));
  if (gk_getwctimer(vault->timer_3) > 0)
    printf("   timer_3: %7.2lfs\n", gk_getwctimer(vault->timer_3));
  if (gk_getwctimer(vault->timer_4) > 0)
    printf("   timer_4: %7.2lfs\n", gk_getwctimer(vault->timer_4));
  if (gk_getwctimer(vault->timer_5) > 0)
    printf("   timer_5: %7.2lfs\n", gk_getwctimer(vault->timer_5));
  if (gk_getwctimer(vault->timer_6) > 0)
    printf("   timer_6: %7.2lfs\n", gk_getwctimer(vault->timer_6));
  printf("----------------------------\n");
  printf("       i/o: %7.2lfs\n", gk_getwctimer(vault->timer_io));
  printf("   tcsetup: %7.2lfs\n", gk_getwctimer(vault->timer_tcsetup));
  printf("  esupport: %7.2lfs\n", gk_getwctimer(vault->timer_esupport));
  printf("   ktsetup: %7.2lfs\n", gk_getwctimer(vault->timer_ktsetup));
  printf(" ktpeeling: %7.2lfs\n", gk_getwctimer(vault->timer_ktpeeling));
  printf("    wclock: %7.2lfs\n", gk_getwctimer(vault->timer_global));
  printf("----------------------------\n");

  printf("RUNRESULT: %s %s %"PRId64" %.2lf %.2lf %.2lf %.2lf %.2lf\n",
      params->infile,
      kttypenames[params->kttype],
      ntriangles,
      gk_getwctimer(vault->timer_tcsetup),
      gk_getwctimer(vault->timer_esupport),
      gk_getwctimer(vault->timer_ktsetup),
      gk_getwctimer(vault->timer_ktpeeling),
      gk_getwctimer(vault->timer_global)
      );

  if(params->stream) {
    printf("Now Streaming....\n\n\n\n");
    char* newEdgesFile = strdup(params->infile);
    char* outputLocation = strdup(params->infile);
    char *ptr = strrchr(newEdgesFile, '/');
    strcpy(ptr, "/newEdges.edges");
    ptr = strrchr(outputLocation, '/');
    strcpy(ptr, "/incremental/");
    // strcpy(ptr, "/");
    printf("Edges file location : %s,     Output Location : %s",
            newEdgesFile, outputLocation);
    stream(params, vault, newEdgesFile, outputLocation);
  }
  else {
    char *ptr = strrchr(params->infile, '/');
    strcpy(ptr, "/baselineTimings.txt");
    FILE *fpout;
    fpout = gk_fopen(params->infile, "a", "fpout");
    fprintf(fpout, "%.6lf\n", gk_getwctimer(vault->timer_global));
    gk_fclose(fpout);
  }

  //cleanMem(vault, params);
}
