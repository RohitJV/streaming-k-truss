#include "kt.h"

#define INF 1000000000

void clearTimers(vault_t *vault) {
  gk_clearwctimer(vault->timer_1);
  gk_clearwctimer(vault->timer_2);
  gk_clearwctimer(vault->timer_3);
  gk_clearwctimer(vault->timer_4);
  gk_clearwctimer(vault->timer_5);
  gk_clearwctimer(vault->timer_6);
  gk_clearwctimer(vault->timer_io);
  gk_clearwctimer(vault->timer_global);
  gk_clearwctimer(vault->timer_esupport);
  gk_clearwctimer(vault->timer_tcsetup);
  gk_clearwctimer(vault->timer_ktsetup);
  gk_clearwctimer(vault->timer_ktpeeling);
}

/*
Write the k-truss number of the entire graph into a file
*/
void writeKTToFile(vault_t *vault,
                   gk_graph_t *graph,
                   char* outputFile) {
  FILE *fpout;
  int32_t nvtxs = graph->nvtxs;
  ssize_t *xadj = graph->xadj;
  int32_t *adjncy = graph->adjncy;
  int32_t *kt = vault->kt;
  ssize_t *endIdx = vault->endIdx;
  int32_t vi, ei;

  fpout = gk_fopen(outputFile, "w", "fpout");
  printf("\nWriting the KT to file : %s\n", outputFile);
  for(vi=0;vi<nvtxs;vi++) {
    for(ei=xadj[vi];ei<endIdx[vi];ei++) {
      fprintf(fpout, "%7d %7d %4d\n",
              vault->iperm[vi]+1, vault->iperm[adjncy[ei]]+1, kt[ei]+2);
    }
  }
  gk_fclose(fpout);
}

/*
Modify the graph for easier traversal
*/
gk_graph_t* createModifiedGraph(vault_t *vault,
                                int buffer) {
  gk_graph_t *original_ugraph = vault->ugraph;
  int32_t nvtxs  = original_ugraph->nvtxs;
  ssize_t *uxadj = original_ugraph->xadj, *xadj;
  int32_t *uadjncy = original_ugraph->adjncy, *adjncy, *kt;

  gk_graph_t *graph = gk_graph_Create();
  graph->nvtxs  = nvtxs;
  graph->xadj   = xadj = gk_zmalloc(nvtxs+1, "createModifiedGraph: modified xadj");
  int32_t modifiedAdjSpace = (uxadj[nvtxs]*2) + buffer*(nvtxs+1);
  graph->adjncy = adjncy = gk_i32malloc(modifiedAdjSpace, "createModifiedGraph: adjncy");
  vault->kt = kt = gk_i32malloc(modifiedAdjSpace, "createModifiedGraph: kt");

  /* Count the length of the adjacency lists, and add buffer space */
  ssize_t ei, vi, vj;
  memset(xadj, 0, (nvtxs+1)*sizeof(ssize_t));
  for(vi=0; vi<nvtxs; vi++) {
    for(ei=uxadj[vi]; ei<uxadj[vi+1]; ei++) {
      vj = uadjncy[ei];
      xadj[vi]++;
      xadj[vj]++;
    }
    xadj[vi] = xadj[vi] + (vi ? xadj[vi-1] : 0) + buffer;
  }

  ssize_t *endIdx = vault->endIdx = gk_zmalloc(nvtxs+1, "createModifiedGraph: start indices of modified xadj");
  /* Fill in the new adjacency list */
  ssize_t *revIdx = vault->revIdx = gk_zmalloc(modifiedAdjSpace, "createModifiedGraph: rev indices of modified xadj");
  for(vi=nvtxs; vi>0; vi--) {
    xadj[vi] = xadj[vi-1];
    endIdx[vi] = xadj[vi];
  }
  xadj[0]=0;
  endIdx[0] = 0;
  for(vi=0; vi<nvtxs; vi++) {
    for(ei=uxadj[vi]; ei<uxadj[vi+1]; ei++) {
      vj = uadjncy[ei];
      /* populate adjancency list */
      adjncy[endIdx[vi]] = vj;
      adjncy[endIdx[vj]] = vi;
      /* Create reverse index for faster lookup later */
      revIdx[endIdx[vi]] = endIdx[vj];
      revIdx[endIdx[vj]] = endIdx[vi];
      /* assert statements - checkpoint */
      if(vault->ktedges[ei].vi != gk_min(vault->iperm[vi], vault->iperm[vj]))
        printf("Something went wrong\n");
      if(vault->ktedges[ei].vj != gk_max(vault->iperm[vi], vault->iperm[vj]))
        printf("Something went wrong.\n");
      /* populate k-truss values */
      kt[endIdx[vi]] = vault->ktedges[ei].k-2;
      kt[endIdx[vj]] = vault->ktedges[ei].k-2;
      endIdx[vi]++; endIdx[vj]++;
    }
  }
  return graph;
}

int32_t addNewEdge(vault_t *vault,
                gk_graph_t *graph,
                edge_t *inserted_edge) {

  int32_t *adjncy=graph->adjncy, *kt = vault->kt;
  ssize_t *revIdx = vault->revIdx;
  int32_t vi=inserted_edge->vi, vj=inserted_edge->vj, retEdge;
  int32_t starti = graph->xadj[vi], endi = vault->endIdx[vi]++;
  int32_t startj = graph->xadj[vj], endj = vault->endIdx[vj]++;

  while(--endi >= starti) {
    if(adjncy[endi]<vj) {
      adjncy[++endi] = vj;
      kt[endi] = vault->ktmax+10;
      retEdge = endi;
      break;
    }
    adjncy[endi+1]=adjncy[endi];
    kt[endi+1]=kt[endi];
    revIdx[revIdx[endi]] = endi+1;
    revIdx[endi+1]=revIdx[endi];
  }
  if(endi<starti) {
    adjncy[++endi] = vj;
    kt[endi] = vault->ktmax+10;
    retEdge = endi;
  }

  while(--endj >= startj) {
    if(adjncy[endj]<vi) {
      adjncy[++endj] = vi;
      kt[endj] = vault->ktmax+10;
      break;
    }
    adjncy[endj+1]=adjncy[endj];
    kt[endj+1]=kt[endj];
    revIdx[revIdx[endj]] = endj+1;
    revIdx[endj+1]=revIdx[endj];
  }
  if(endj<startj) {
    adjncy[++endj] = vi;
    kt[endj] = vault->ktmax+10;
  }

  /* Update revIdx */
  revIdx[endj] = endi;
  revIdx[endi] = endj;

  endi = vault->endIdx[vi];
  endj = vault->endIdx[vj];

  /* return one edge index for later use */
  return retEdge;
}

/*
Returns the number of triangles on an edge 'e', such that the truss number of
those triangles is atleast K(e).
*/
int32_t MTD(vault_t *vault,
            gk_graph_t *graph,
            int32_t *mtd,
            int32_t e) {

  if(mtd[e]!=-1)
    return mtd[e];
  mtd[e]=0;

  int32_t *adjncy=graph->adjncy, *kt = vault->kt;
  ssize_t *revIdx = vault->revIdx;
  int32_t vi=adjncy[e], vj=adjncy[revIdx[e]], ei, ej, k=kt[e];
  int32_t starti = graph->xadj[vi], endi = vault->endIdx[vi];
  int32_t startj = graph->xadj[vj], endj = vault->endIdx[vj];

  while(starti<endi && startj<endj) {
    if(adjncy[starti] == adjncy[startj]) {
      ei = starti;
      ej = startj;
      if(gk_min(kt[ei], kt[ej]) >= k)
        mtd[e]++;
      starti++; startj++;
    }
    else if(adjncy[starti] < adjncy[startj])
      starti++;
    else
      startj++;
  }

  return (mtd[revIdx[e]] = mtd[e]);
}

/*
Returns the number of triangles on an edge 'e', such that
1) The truss number of the triangle is greater that K(e), or
2) The truss number of the triangle is equal to K(e) and the MTD(e) is greater
   than K(e).
*/
int32_t PTD(vault_t *vault,
            gk_graph_t *graph,
            int32_t *ptd,
            int32_t *mtd,
            int32_t e) {

  if(ptd[e]!=-1)
    return ptd[e];
  ptd[e]=0;

  int32_t *adjncy=graph->adjncy, *kt = vault->kt;
  ssize_t *revIdx = vault->revIdx;
  int32_t vi=adjncy[e], vj=adjncy[revIdx[e]], ei, ej, k=kt[e];
  int32_t starti = graph->xadj[vi], endi = vault->endIdx[vi];
  int32_t startj = graph->xadj[vj], endj = vault->endIdx[vj];

  while(starti<endi && startj<endj) {
    if(adjncy[starti] == adjncy[startj]) {
      ei = starti;
      ej = startj;
      if(gk_min(kt[ei], kt[ej]) > k)
        ptd[e]++;
      else if(gk_min(kt[ei], kt[ej]) == k) {
        if(kt[ei]!=kt[ej]) {
          if( (kt[ei]<kt[ej] && MTD(vault, graph, mtd, ei)>k) || (kt[ej]<kt[ei] && MTD(vault, graph, mtd, ej)>k) )
            ptd[e]++;
        }
        else {
          if(MTD(vault, graph, mtd, ei)>k && MTD(vault, graph, mtd, ej)>k)
            ptd[e]++;
        }
      }
      starti++; startj++;
    }
    else if(adjncy[starti] < adjncy[startj])
      starti++;
    else
      startj++;
  }

  return (ptd[revIdx[e]] = ptd[e]);
}

void stream(params_t *params,
            vault_t *vault,
            char* edgesFile,
            char* outputLocation) {

  int32_t buffer = 1000; /* How many edges (max) can it handle */
  /* Initialize timers */
  double totalRuntime = 0.0, addEdgeTime = 0.0, selectRootEdgesTime = 0.0, traversalTime = 0.0;

  /* Pre-processing done here */
  gk_graph_t *modifiedGraph = createModifiedGraph(vault, buffer);
  ssize_t nvtxs = vault->graph->nvtxs;
  edge_t* inserted_edge = gk_malloc(sizeof(edge_t), "stream: inserted edge");
  int32_t *rootEdges = gk_i32malloc(vault->graph->nvtxs * 2, "stream: Root edges to begin the search space from");
  int32_t *visited = gk_i32malloc(modifiedGraph->xadj[nvtxs], "stream: visited");
  int32_t *evicted = gk_i32malloc(modifiedGraph->xadj[nvtxs], "stream: evicted");
  int32_t *ptd = gk_i32malloc(modifiedGraph->xadj[nvtxs], "stream: pure truss degree");
  int32_t *mtd = gk_i32malloc(modifiedGraph->xadj[nvtxs], "stream: maximum truss degree");
  int32_t *td = gk_i32malloc(modifiedGraph->xadj[nvtxs], "stream: prospective truss degree");

  /* Read edges file*/
  FILE *fpin=NULL;
  size_t lnlen;
  char *line=NULL;
  int number_of_edges;
  fpin = gk_fopen(edgesFile, "r", "stream: fpin");
  do {
    if (gk_getline(&line, &lnlen, fpin) <= 0)
      gk_errexit(SIGERR, "Premature end of input file: file\n");
  } while (line[0] == '%');
  sscanf(line, "%d", &number_of_edges);

  int i=0;
  buffer = number_of_edges;

  for(int i=0; i<number_of_edges; i++) {
    do {
      if (gk_getline(&line, &lnlen, fpin) == -1)
        gk_errexit(SIGERR, "Pregraphure end of input file: file while reading row %d\n", i);
    } while (line[0] == '%');
    int vtx1, vtx2;
    sscanf(line, "%d %d", &vtx1, &vtx2);

    clearTimers(vault);
    /* Initialize memory */
    memset(visited, 0, modifiedGraph->xadj[nvtxs] * sizeof(int32_t));
    memset(evicted, 0, modifiedGraph->xadj[nvtxs] * sizeof(int32_t));
    memset(mtd, -1, modifiedGraph->xadj[nvtxs] * sizeof(int32_t));
    memset(ptd, -1, modifiedGraph->xadj[nvtxs] * sizeof(int32_t));
    memset(td, 0, modifiedGraph->xadj[nvtxs] * sizeof(int32_t));

    gk_startwctimer(vault->timer_global);
    printf("Inserted Edge %d: (%d, %d)\n", i, vtx1-1, vtx2-1);
    inserted_edge->vi = vault->perm[vtx1-1];
    inserted_edge->vj = vault->perm[vtx2-1];

    /* Add the new edge to the modifiedGraph */
    gk_startwctimer(vault->timer_1);
    int32_t refEdge = addNewEdge(vault, modifiedGraph, inserted_edge);
    gk_stopwctimer(vault->timer_1);

    /* Select rootedges to begin DFS from */
    gk_startwctimer(vault->timer_2);
    int32_t *trussWiseSupCount = gk_i32malloc(vault->ktmax+2, "stream: Truss-wise support count");
    int32_t rootEdgesCount = 0;
    selectRootEdgesForDFS(vault, modifiedGraph, inserted_edge, trussWiseSupCount, rootEdges, &rootEdgesCount);
    gk_stopwctimer(vault->timer_2);

    /* Traversal algorithm */
    gk_startwctimer(vault->timer_3);
    int evictionTime = 1;
    traversalAlgorithm(vault, modifiedGraph, refEdge, visited, evicted,
                       trussWiseSupCount, rootEdges, rootEdgesCount,
                       ptd, mtd, td, &evictionTime);
    gk_stopwctimer(vault->timer_3);
    gk_stopwctimer(vault->timer_global);

    /* Write to output file */
#ifdef TEST_CORRECTNESS
    char* outputFile = strdup(outputLocation);
    char outputNum[10];
    sprintf(outputNum, "%d", i+1);
    strcat(outputFile, outputNum); strcat(outputFile, ".out");
    writeKTToFile(vault, modifiedGraph, outputFile);
#endif

    gk_free((void **)&trussWiseSupCount, LTERM);
    /* Calculate runtimes of various components */
    totalRuntime = totalRuntime + gk_getwctimer(vault->timer_global);
    addEdgeTime = addEdgeTime + gk_getwctimer(vault->timer_1);
    selectRootEdgesTime = selectRootEdgesTime + gk_getwctimer(vault->timer_2);
    traversalTime = traversalTime + gk_getwctimer(vault->timer_3);

    char *ptr = strrchr(params->infile, '/');
    strcpy(ptr, "/incrementalTimings.txt");
    FILE *fpout;
    fpout = gk_fopen(params->infile, "a", "fpout");
    fprintf(fpout, "%.6lf\n", gk_getwctimer(vault->timer_global));
    gk_fclose(fpout);
  }

  printf("Total addEdgeTime Runtime : %.2lf\n", addEdgeTime);
  printf("Total selectRootEdgesTime Runtime : %.2lf\n", selectRootEdgesTime);
  printf("Total Incremental traversalTime : %.2lf\n", traversalTime);
  printf("Total Incremental Runtime : %.2lf\n", totalRuntime);

  /* Cleanup */
  gk_free((void **)&inserted_edge, LTERM);
  gk_free((void **)&rootEdges, LTERM);
  gk_free((void **)&visited, LTERM);
  gk_free((void **)&evicted, LTERM);
  gk_free((void **)&mtd, LTERM);
  gk_free((void **)&ptd, LTERM);
  gk_free((void **)&td, LTERM);
  gk_free((void **)&vault->kt, LTERM);
  gk_free((void **)&vault->endIdx, LTERM);
  gk_free((void **)&vault->revIdx, LTERM);
  gk_free((void **)&modifiedGraph->xadj, LTERM);
  gk_free((void **)&modifiedGraph->adjncy, LTERM);
  gk_free((void **)&modifiedGraph, LTERM);
}
