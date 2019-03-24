#include "kt.h"
#include <stdarg.h>

#define STATIC_ADJ 0
#define STREAM_ADJ 1

#define FIND_LOCAL_KTMAX 10
#define SELECT_ROOT_EDGES 11

typedef struct {
  int32_t *rootEdgeLocations;
  eedge_t **rootEdges;
  int8_t ktmax;
} selectRootEdges_struct;

static void clearTimers(vault_t *vault) {
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

static ggraph_t* createModifiedGraph(vault_t *vault) {
  gk_graph_t *original_ugraph = vault->ugraph;
  ssize_t *uxadj = original_ugraph->xadj;
  int32_t *uadjncy = original_ugraph->adjncy;

  /* Initialize the new graph */
  ggraph_t *modifiedGraph = gk_malloc(sizeof(ggraph_t), "createModifiedGraph: graph for streaming");
  int32_t nvtxs = modifiedGraph->nvtxs = original_ugraph->nvtxs;
  ssize_t *xadj = modifiedGraph->xadj = gk_zsmalloc(nvtxs+1, 0, "createModifiedGraph: modified xadj");
  eedge_t *adjncy = modifiedGraph->edges = gk_malloc(sizeof(eedge_t)*(uxadj[nvtxs]*2+1), "createModifiedGraph: adjacency list");

  /* Count the length of the adjacency lists */
  ssize_t e, vi, vj;
  for(vi=0; vi<nvtxs; vi++) {
    for(e=uxadj[vi]; e<uxadj[vi+1]; e++) {
      vj = uadjncy[e];
      xadj[vi]++;
      xadj[vj]++;
    }
    xadj[vi] = xadj[vi] + (vi ? xadj[vi-1] : 0);
  }

  /* Shift the xadj */
  ssize_t *endIdx = gk_zmalloc(nvtxs+1, "createModifiedGraph: temporary ds to populate the graph");
  for(vi=nvtxs; vi>0; vi--) {
    xadj[vi] = xadj[vi-1];
    endIdx[vi] = xadj[vi];
  }
  endIdx[0] = xadj[0] = 0;

  /* Populate the graph */
  for(vi=0; vi<nvtxs; vi++) {
    for(e=uxadj[vi]; e<uxadj[vi+1]; e++) {
      vj = uadjncy[e];
      int32_t ei=endIdx[vi]++, ej=endIdx[vj]++;
      adjncy[ei].v = vj; adjncy[ej].v = vi;
      adjncy[ei].revPtr = &adjncy[ej]; adjncy[ej].revPtr = &adjncy[ei];
      adjncy[ei].k = adjncy[ej].k = vault->ktedges[e].k-2;
      adjncy[ei].mtd = adjncy[ei].ptd = 1; adjncy[ej].mtd = adjncy[ej].ptd = -1;
      adjncy[ei].td = adjncy[ei].visited = adjncy[ej].td = adjncy[ej].visited = 0;
    }
  }

  /* Free resources */
  gk_free((void **)&endIdx, LTERM);
  gk_free((void **)&vault->graph->xadj, LTERM);
  gk_free((void **)&vault->graph->adjncy, LTERM);
  gk_free((void **)&vault->ugraph->xadj, LTERM);
  gk_free((void **)&vault->ugraph->adjncy, LTERM);
  gk_free((void **)&vault->ugraph->iadjwgt, LTERM);
  gk_free((void **)&vault->ugraph, LTERM);
  gk_free((void **)&vault->lgraph->xadj, LTERM);
  gk_free((void **)&vault->lgraph->adjncy, LTERM);
  gk_free((void **)&vault->lgraph, LTERM);
  gk_free((void **)&vault->ktedges, LTERM);

  return modifiedGraph;
}

static int32_t* readNewEdges(char* edgesFile, ggraph_t *modifiedGraph, vault_t *vault) {
  size_t lnlen;
  char *line=NULL;
  int32_t number_of_edges;
  FILE *fpin = gk_fopen(edgesFile, "r", "stream: fpin");
  do {
    if (gk_getline(&line, &lnlen, fpin) <= 0)
      gk_errexit(SIGERR, "Premature end of input file: file\n");
  } while (line[0] == '%');
  sscanf(line, "%d", &number_of_edges);

  int32_t *newEdges = gk_i32malloc(1 + number_of_edges*2, "readNewEdges: vertices of the edges");
  newEdges[0] = number_of_edges;
  for(int32_t i=0; i<number_of_edges; i++) {
    do {
      if (gk_getline(&line, &lnlen, fpin) == -1)
        gk_errexit(SIGERR, "Pregraphure end of input file: file while reading row %d\n", i);
    } while (line[0] == '%');
    int32_t vtx1, vtx2;
    sscanf(line, "%d %d", &vtx1, &vtx2);
    /* Place the new edge vertices as vi followed by vj, with vi<vj */
    newEdges[i*2+1] = gk_min(vault->perm[vtx1-1], vault->perm[vtx2-1]);
    newEdges[i*2+2] = gk_max(vault->perm[vtx1-1], vault->perm[vtx2-1]);
  }

  return newEdges;
}

static void insertToHmapChain(edgeblock_t **head, edgeblock_t *edgeBlock) {
  edgeblock_t *prev=NULL, *cur=*head;
  int32_t v = edgeBlock->edge.v;
  while(cur!=NULL && cur->edge.v < v) {
    prev=cur;
    cur=cur->next;
  }
  if(prev!=NULL) {
    edgeBlock->next = prev->next;
    prev->next=edgeBlock;
  }
  else {
    edgeBlock->next = *head;
    *head=edgeBlock;
  }
}

static void addNewEdge(edgeblock_t **hmap, int32_t vi, int32_t vj) {
  int32_t hval;

  /* Find the correct location in hmap to place the new edge */
  for(hval = vi&(HMAP_SIZE-1); hmap[hval]!=NULL && hmap[hval]->src!=vi; hval=(hval+1)&(HMAP_SIZE-1));
  /* Place the new edge into chain at the hmap location */
  edgeblock_t *edgeBlock1 = gk_malloc(sizeof(edgeblock_t), "addNewEdge: create a new edgeblock");
  edgeBlock1->src = vi; edgeBlock1->edge.v=vj; edgeBlock1->next=NULL;
  if(hmap[hval]==NULL)
    hmap[hval]=edgeBlock1;
  else
    insertToHmapChain(&hmap[hval], edgeBlock1);

/* Repeat above for the reverse edge */
  /* Find the correct location in hmap to place the new edge */
  for(hval = vj&(HMAP_SIZE-1); hmap[hval]!=NULL && hmap[hval]->src!=vj; hval=(hval+1)&(HMAP_SIZE-1));
  /* Place the new edge into chain at the hmap location */
  edgeblock_t *edgeBlock2 = gk_malloc(sizeof(edgeblock_t), "addNewEdge: create a new edgeblock");
  edgeBlock2->src = vj; edgeBlock2->edge.v=vi; edgeBlock2->next=NULL;
  if(hmap[hval]==NULL)
    hmap[hval]=edgeBlock2;
  else
    insertToHmapChain(&hmap[hval], edgeBlock2);

  /* Cross reference edgeBlock1 and edgeBlock2; populate other details */
  edgeBlock1->edge.revPtr = &(edgeBlock2->edge); edgeBlock2->edge.revPtr = &(edgeBlock1->edge);
  edgeBlock1->edge.k = edgeBlock1->edge.td = edgeBlock1->edge.visited = 0;
  edgeBlock2->edge.k = edgeBlock2->edge.td = edgeBlock2->edge.visited = 0;
  edgeBlock1->edge.mtd = edgeBlock1->edge.ptd = -1;
  edgeBlock2->edge.mtd = edgeBlock2->edge.ptd = -1;
}

static edgeblock_t* findMapEntry(edgeblock_t **hmap, int32_t vtx) {
  int32_t hval;
  for(hval = vtx&(HMAP_SIZE-1); hmap[hval]!=NULL && hmap[hval]->src!=vtx; hval=(hval+1)&(HMAP_SIZE-1));
  return hmap[hval];
}

static void intersectLists(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap,
                           int32_t vi, int32_t vj, void* result, int8_t optype) {

  int32_t starti = modifiedGraph->xadj[vi], startj = modifiedGraph->xadj[vj];
  int32_t endi = modifiedGraph->xadj[vi+1], endj = modifiedGraph->xadj[vj+1];
  eedge_t *adjncy = modifiedGraph->edges;
  edgeblock_t *ptri=findMapEntry(hmap, vi), *ptrj=findMapEntry(hmap, vj);
  eedge_t *ei, *ej;

  while((starti<endi || ptri!=NULL) && (startj<endj || ptrj!=NULL)) { /* This guarantees that we have reached the end of the actual adjancency lists */
    /* Get the actual current entry from vi's adjacency list */
    int8_t flagi, flagj, flag = 3;
    if(ptri!=NULL && starti<endi) {
      if(ptri->edge.v < adjncy[starti].v) {
        ei = &(ptri->edge);
        flagi = STREAM_ADJ;
      }
      else {
        ei = &adjncy[starti];
        flagi = STATIC_ADJ;
      }
    }
    else if(ptri!=NULL) {
      ei = &(ptri->edge);
      flagi = STREAM_ADJ;
    }
    else {
      ei = &adjncy[starti];
      flagi = STATIC_ADJ;
    }
    /* Repeat for vj */
    if(ptrj!=NULL && startj<endj) {
      if(ptrj->edge.v < adjncy[startj].v) {
        ej = &(ptrj->edge);
        flagj = STREAM_ADJ;
      }
      else {
        ej = &adjncy[startj];
        flagj = STATIC_ADJ;
      }
    }
    else if(ptrj!=NULL) {
      ej = &(ptrj->edge);
      flagj = STREAM_ADJ;
    }
    else {
      ej = &adjncy[startj];
      flagj = STATIC_ADJ;
    }

    /* Perform operations on the 2 edges */
    if(optype == FIND_LOCAL_KTMAX) {
      int32_t* trussWiseSupCount = (int32_t *)result;
      int32_t minK = gk_min(ei->k, ej->k);
      trussWiseSupCount[minK]++;
    }
    else if(optype == SELECT_ROOT_EDGES) {
      selectRootEdges_struct *selectRootEdgesStruct = (selectRootEdges_struct *)result;
      int32_t *rootEdgeLocations = selectRootEdgesStruct->rootEdgeLocations;
      eedge_t **rootEdges = selectRootEdgesStruct->rootEdges;
      int8_t local_ktmax = selectRootEdgesStruct->ktmax;
      int32_t minK = gk_min(ei->k, ej->k);
      if(minK <= local_ktmax)
        rootEdges[--rootEdgeLocations[minK]] = (ei->k < ej->k ? ei : ej);
    }

    /* Increment edge pointers appropriately. */
    int32_t vali = ei->v, valj = ej->v;
    if(vali < valj) /* Increment ei */
      flag = 2;
    else if(vali > valj) /* Increment ej */
      flag = 1;
    if(flag & 2) { /* Increment ei */
      if(flagi == STATIC_ADJ)
        starti++;
      else
        ptri=ptri->next;
    }
    if(flag & 1) { /* Increment ej */
      if(flagj == STATIC_ADJ)
        startj++;
      else
        ptrj=ptrj->next;
    }
  }
}

static void selectRootEdges(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap, int32_t vi, int32_t vj) {
  /* Populate the support-count array (trussWiseSupCount[]) */
  int32_t *trussWiseSupCount = gk_i32smalloc(vault->ktmax+1, 0, "selectRootEdges: truss-wise support count");
  int32_t *trussWiseSupCountSuffixSum = gk_i32smalloc(vault->ktmax+1, 0, "selectRootEdges: truss-wise support count with suffix sum");
  intersectLists(vault, modifiedGraph, hmap, vi, vj, (void *)trussWiseSupCount, FIND_LOCAL_KTMAX);

  /* Calculate the local_ktmax */
  trussWiseSupCountSuffixSum[vault->ktmax] = trussWiseSupCount[vault->ktmax];
  for(int8_t i=vault->ktmax-1; i>=0; i--)
    trussWiseSupCountSuffixSum[i] = trussWiseSupCount[i] + trussWiseSupCountSuffixSum[i+1];
  int8_t local_ktmax;
  for(local_ktmax=vault->ktmax; local_ktmax>=0; local_ktmax--) {
    if(trussWiseSupCountSuffixSum[local_ktmax] <= local_ktmax)
      continue;
    break;
  }

  /* Populate the root edges with 0 <= K(.) <= local_ktmax */
  selectRootEdges_struct selectRootEdgesStruct;
  selectRootEdgesStruct.rootEdgeLocations = trussWiseSupCountSuffixSum;
  selectRootEdgesStruct.rootEdges = gk_malloc(modifiedGraph->nvtxs * 2 * sizeof(eedge_t*), "selectRootEdges: pointers to root edges");
  selectRootEdgesStruct.ktmax = local_ktmax;
  intersectLists(vault, modifiedGraph, hmap, vi, vj, (void *)&selectRootEdgesStruct, SELECT_ROOT_EDGES);

  /* Free Resources */
  gk_free((void **)&trussWiseSupCount, LTERM);
  gk_free((void **)&trussWiseSupCountSuffixSum, LTERM);
}

static void printAdjList(ggraph_t *modifiedGraph, edgeblock_t **hmap, int32_t v) {
  printf("\nAdjancency list of the vertex : %d\n", v);
  int32_t start = modifiedGraph->xadj[v], end = modifiedGraph->xadj[v+1];
  eedge_t *adjncy = modifiedGraph->edges;
  edgeblock_t *ptr=findMapEntry(hmap, v);

  while(start < end) {
    printf("%d(o), ", adjncy[start].v);
    start++;
  }
  while(ptr!=NULL) {
    printf("%d(+), ", ptr->edge.v);
    ptr=ptr->next;
  }
}

void stream_v2(params_t *params, vault_t *vault, char* edgesFile, char* outputLocation) {
  double totalRuntime = 0.0, addEdgeTime = 0.0, selectRootEdgesTime = 0.0, traversalTime = 0.0;

  /* Pre-processing done here */
  ggraph_t *modifiedGraph = createModifiedGraph(vault); /* Create the new graph */
  edgeblock_t *hmap[HMAP_SIZE];
  for(int32_t i=0; i<HMAP_SIZE; i++)
    hmap[i]=NULL;

  /* Read new edges from file */
  int32_t *newEdges = readNewEdges(edgesFile, modifiedGraph, vault);
  int32_t idx, number_of_edges = newEdges[0];
  for(idx=0; idx<number_of_edges; idx++) {
    /* Process one edge at a time */
    clearTimers(vault);
    gk_startwctimer(vault->timer_global);
    int32_t vi=newEdges[idx*2+1], vj=newEdges[idx*2+2];

    printf("\nNew Edge : (%d, %d)\n", vi, vj);
    gk_startwctimer(vault->timer_1);
    /* 1. Add a new edge to the hashmap */
    addNewEdge(hmap, vi, vj);
    // printAdjList(modifiedGraph, hmap, vi);
    // printAdjList(modifiedGraph, hmap, vj);
    gk_stopwctimer(vault->timer_1);

    /* 2. Find root edges */
    gk_startwctimer(vault->timer_2);
    selectRootEdges(vault, modifiedGraph, hmap, vi, vj);
    gk_stopwctimer(vault->timer_2);

    gk_stopwctimer(vault->timer_global);
    /* Output the runtimes to file */
    addEdgeTime = addEdgeTime + gk_getwctimer(vault->timer_1);
    selectRootEdgesTime = selectRootEdgesTime + gk_getwctimer(vault->timer_2);
    traversalTime = traversalTime + gk_getwctimer(vault->timer_3);
    totalRuntime = totalRuntime + gk_getwctimer(vault->timer_global);
    char *ptr = strrchr(params->infile, '/');
    strcpy(ptr, "/incrementalTimings.txt");
    FILE *fpout;
    fpout = gk_fopen(params->infile, "a", "fpout");
    fprintf(fpout, "%.6lf\n", gk_getwctimer(vault->timer_global));
    gk_fclose(fpout);
  }

  printf("\n\n\nTotal addEdgeTime Runtime : %.8lf\n", addEdgeTime);
  printf("Total selectRootEdgesTime Runtime : %.8lf\n", selectRootEdgesTime);
  printf("Total Incremental traversalTime : %.8lf\n", traversalTime);
  printf("Total Incremental Runtime : %.8lf\n", totalRuntime);

  /* Cleanup */
  gk_free((void **)&newEdges, LTERM);
  gk_free((void **)&vault->kt, LTERM);
  gk_free((void **)&vault->endIdx, LTERM);
  gk_free((void **)&vault->revIdx, LTERM);
}
