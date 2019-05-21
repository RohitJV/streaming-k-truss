#include "kt.h"
#include <stdarg.h>

#define STATIC_ADJ 0
#define STREAM_ADJ 1

#define SELECT_ROOT_EDGES 10
#define EXPLORE_SEARCHSPACE 11
#define EVICT_EDGES 12
#define UPDATE_MTD 13
#define INSERT_EDGE_MTD 14

#define QUEUE_SIZE 10
#define SEARCHSPACE_QUEUE_SIZE 1000

typedef struct {
  int32_t front, rear, size;
  eedge_t **edgelist;
} edge_q;

typedef struct {
  edge_q *queue;
  int32_t ktmax;
} rootEdges_struct;

typedef struct {
  edge_q *queue;
  eedge_t *e;
  int32_t k;
} searchSpace_struct;

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

static void writeKTToFile(vault_t *vault, ggraph_t *graph, edgeblock_t **hmap, char* outputFile) {
  FILE *fpout;
  int32_t nvtxs = graph->nvtxs;
  ssize_t *xadj = graph->xadj;
  eedge_t *adjncy = graph->edges;
  ssize_t *endIdx = vault->endIdx;
  int32_t vi, ei;

  fpout = gk_fopen(outputFile, "w", "fpout");
  printf("Writing the KT to file : %s\n", outputFile);
  /* Print the truss numbers of the static graph */
  for(vi=0;vi<nvtxs;vi++) {
    for(ei=xadj[vi];ei<xadj[vi+1];ei++) {
      fprintf(fpout, "%7d %7d %4d\n",
              vault->iperm[vi]+1, vault->iperm[adjncy[ei].v]+1, adjncy[ei].k+2);
    }
  }
  /* Print the truss numbers of the inserted edges */
  for(int32_t i=0; i<HMAP_SIZE; i++) {
    edgeblock_t *cur=hmap[i];
    while(cur!=NULL) {
	    fprintf(fpout, "%7d %7d %4d\n",
              vault->iperm[cur->edge.v]+1, vault->iperm[(cur->edge.revPtr)->v]+1, cur->edge.k+2);
      cur=cur->next;
    }
  }
  gk_fclose(fpout);
}

static edgeblock_t* findMapEntry(edgeblock_t **hmap, int32_t vtx) {
  int32_t hval;
  for(hval = vtx&(HMAP_SIZE-1); hmap[hval]!=NULL && hmap[hval]->src!=vtx; hval=(hval+1)&(HMAP_SIZE-1));
  return hmap[hval];
}

static int64_t intersectLists(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap,
                           int32_t vi, int32_t vj, void* result, int8_t optype) {

  int32_t starti = modifiedGraph->xadj[vi], startj = modifiedGraph->xadj[vj];
  int32_t endi = modifiedGraph->xadj[vi+1], endj = modifiedGraph->xadj[vj+1];
  eedge_t *adjncy = modifiedGraph->edges;
  edgeblock_t *ptri=findMapEntry(hmap, vi), *ptrj=findMapEntry(hmap, vj);
  eedge_t *ei, *ej;
  int64_t ntriangles = 0;
  int32_t count = 0;

  while((starti<endi || ptri!=NULL) && (startj<endj || ptrj!=NULL)) { /* This guarantees that we have reached the end of the actual adjancency lists */
    ntriangles++;
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

    /* Perform operations on the 2 edges, if they share a common vertex*/
    if(ei->v == ej->v) {
      /* Pick the edge with the lower K(.) */
      int8_t minK;
      eedge_t *minE, *theOtherE=NULL;
      if(ei->k < ej->k) {
        minK=ei->k;
        minE=ei;
      }
      else if(ei->k == ej->k) {
        minK=ei->k;
        minE=ei;
        theOtherE=ej;
      }
      else {
        minK=ej->k;
        minE=ej;
      }

      if(optype == SELECT_ROOT_EDGES) {
        minE->mtd = minE->mtd + 1;
        minE->revPtr->mtd = minE->revPtr->mtd + 1;
        if(theOtherE!=NULL) {
          theOtherE->mtd = theOtherE->mtd + 1;
          theOtherE->revPtr->mtd = theOtherE->revPtr->mtd + 1;
        }
        if(minE->mtd>minK && (theOtherE==NULL || theOtherE->mtd>minK)) {
          rootEdges_struct *rootEdges = (rootEdges_struct *)result;
          int32_t pos = ++(rootEdges->queue[minK].rear);
          if(pos==rootEdges->queue[minK].size-1) {
            rootEdges->queue[minK].size = rootEdges->queue[minK].size * 2;
            rootEdges->queue[minK].edgelist = gk_realloc(rootEdges->queue[minK].edgelist, sizeof(eedge_t*) * rootEdges->queue[minK].size,
                                                          "selectRootEdges : realloc queue of edges for a value of K(.)");
          }
          rootEdges->queue[minK].edgelist[pos]=minE;
        }
      }
      else if(optype == EXPLORE_SEARCHSPACE) {
        searchSpace_struct *searchSpace = (searchSpace_struct *)result;
        eedge_t *e = searchSpace->e;
        edge_q *searchQueue = searchSpace->queue;
        if(minK>searchSpace->k || (minK==searchSpace->k && minE->mtd>minK && (theOtherE==NULL || theOtherE->mtd>minK))) {
          e->td++;
          e->revPtr->td++;        
          if(minK==searchSpace->k && !minE->visited) {
            minE->visited = 1;
            minE->revPtr->visited = 1;
            searchQueue->rear++;
            count++;
            if(searchQueue->rear == searchQueue->size) {
              searchQueue->size *= 2;
              searchQueue->edgelist = gk_realloc(searchQueue->edgelist,  sizeof(eedge_t*) * searchQueue->size, "updateTrussNumbers : realloc queue of edges in the initial search space");
            }
            searchQueue->edgelist[searchQueue->rear] = minE;
          }
        	if(minK==searchSpace->k && theOtherE!=NULL && !theOtherE->visited) {
        	  theOtherE->visited = 1;
            theOtherE->revPtr->visited = 1;
            searchQueue->rear++;
            count++;
            if(searchQueue->rear == searchQueue->size) {
              searchQueue->size *= 2;
              searchQueue->edgelist = gk_realloc(searchQueue->edgelist,  sizeof(eedge_t*) * searchQueue->size, "updateTrussNumbers : realloc queue of edges in the initial search space");
            }
            searchQueue->edgelist[searchQueue->rear] = theOtherE;
        	}
        }
      }
      else if(optype == EVICT_EDGES) {
        searchSpace_struct *searchSpace = (searchSpace_struct *)result;
        eedge_t *e = searchSpace->e;
        edge_q *evictionQueue = searchSpace->queue;
        /*
        Consider only the triangles with minK = k, and
        --> If one of the edges of the triangle was already evicted, the other edge already had its truss-degree reduced due to this triangle.
            So, skip reducing the truss-degree of the other edge.
        --> Otherwise, push the edge(s) with K(.)=k to the eviction queue, as long as they don't already exist there.
        */
        if(minK==searchSpace->k && minE->visited!=-2 && (theOtherE==NULL || theOtherE->visited!=-2)) {
          //if(minE->visited) { /* Check if this edge was even present in the searchSpace */
            minE->td--;
            minE->revPtr->td--;
            if(minE->visited!=-1 && minE->td<=searchSpace->k) { /* Insert minE to the eviction queue if it is not already there */
              minE->visited = -1;
              minE->revPtr->visited = -1;
              evictionQueue->rear++;
              if(evictionQueue->rear == evictionQueue->size) {
                evictionQueue->size *= 2;
                evictionQueue->edgelist = gk_realloc(evictionQueue->edgelist,  sizeof(eedge_t*) * evictionQueue->size, "updateTrussNumbers : realloc eviction list");
              }
              evictionQueue->edgelist[evictionQueue->rear] = minE;
            }
          //}
          if(theOtherE!=NULL) {
              //&& theOtherE->visited) { /* Check if this edge was even present in the searchSpace */
            theOtherE->td--;
            theOtherE->revPtr->td--;
            if(theOtherE->visited!=-1 && theOtherE->td<=searchSpace->k) { /* Insert theOtherE to the eviction queue if it is not already there */
              theOtherE->visited=-1;
              theOtherE->revPtr->visited=-1;
              evictionQueue->rear++;
              if(evictionQueue->rear == evictionQueue->size) {
                evictionQueue->size *= 2;
                evictionQueue->edgelist = gk_realloc(evictionQueue->edgelist, sizeof(eedge_t*) *  evictionQueue->size, "updateTrussNumbers : realloc eviction list");
              }
              evictionQueue->edgelist[evictionQueue->rear] = theOtherE;
            }
          }
        }
      }
      else if(optype == UPDATE_MTD) {
        eedge_t* e = (eedge_t*)result;
        if(minK>=e->k)
          e->mtd = e->mtd + 1;
      }
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

  return ntriangles;
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
      adjncy[ei].mtd = adjncy[ej].mtd = 0;
      adjncy[ei].td = adjncy[ei].visited = adjncy[ej].td = adjncy[ej].visited = 0;
    }
  }

  /* Populate the mtd values */
  for(int32_t ei=0; ei<xadj[nvtxs]; ei++) {
    if(adjncy[ei].mtd)
      continue;
    vi = adjncy[ei].v;
    vj = (adjncy[ei].revPtr)->v;
    int8_t k = adjncy[ei].k;
    int32_t starti=xadj[vi], endi=endIdx[vi], startj=xadj[vj], endj=endIdx[vj];
    while(starti<endi && startj<endj) {
      if(adjncy[starti].v == adjncy[startj].v) {
        if(gk_min(adjncy[starti].k, adjncy[startj].k) >= k) {
          adjncy[ei].mtd = adjncy[ei].mtd + 1;
          adjncy[ei].revPtr->mtd = adjncy[ei].revPtr->mtd + 1;
        }
        starti++;
        startj++;
      }
      else if(adjncy[starti].v < adjncy[startj].v)
        starti++;
      else
        startj++;
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

static edgeblock_t* addNewEdge(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap, int32_t vi, int32_t vj) {
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
  edgeBlock1->edge.k = edgeBlock2->edge.k = vault->ktmax + 10;
  edgeBlock1->edge.td = edgeBlock1->edge.visited = 0;
  edgeBlock2->edge.td = edgeBlock2->edge.visited = 0;
  edgeBlock1->edge.mtd = 0; edgeBlock2->edge.mtd = 0;

  int64_t ntriangles = intersectLists(vault, modifiedGraph, hmap, vi, vj, NULL, INSERT_EDGE_MTD);

  return edgeBlock1;
}

static rootEdges_struct* selectRootEdges(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap, int32_t vi, int32_t vj) {
  rootEdges_struct *rootEdges = gk_malloc(sizeof(rootEdges_struct), "selectRootEdges : rootEdges_struct");
  rootEdges->ktmax = vault->ktmax;
  rootEdges->queue = gk_malloc(sizeof(edge_q) * (vault->ktmax+1), "selectRootEdges : array of queue of edges for each value of K(.)");
  for(int32_t i=0; i<=vault->ktmax; i++) {
    rootEdges->queue[i].front = rootEdges->queue[i].rear = -1;
    rootEdges->queue[i].size = QUEUE_SIZE;
    rootEdges->queue[i].edgelist = gk_malloc(sizeof(eedge_t*) * QUEUE_SIZE, "selectRootEdges : queue of edges for a value of K(.)");
  }
  int64_t ntriangles = intersectLists(vault, modifiedGraph, hmap, vi, vj, (void *)rootEdges, SELECT_ROOT_EDGES);
  int32_t suffixsum = 0;
  int32_t local_ktmax;
  for(local_ktmax=vault->ktmax; local_ktmax>=0; local_ktmax--) {
    suffixsum = suffixsum + rootEdges->queue[local_ktmax].rear + 1;
    if(suffixsum>local_ktmax)
      break;
    gk_free((void **)&rootEdges->queue[local_ktmax].edgelist, LTERM);
  }
  rootEdges->ktmax = local_ktmax;
  return rootEdges;
}

static void updateMTD(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap, edgeblock_t* inserted_edge) {
  /* Populate the mtd values */
  int32_t ntriangles;
  int32_t nvtxs = modifiedGraph->nvtxs;
  eedge_t *adjncy = modifiedGraph->edges;
  ssize_t *xadj = modifiedGraph->xadj;
  ssize_t e, vi, vj;
  for(int32_t ei=0; ei<xadj[nvtxs]; ei++) {
    vi = adjncy[ei].v;
    vj = (adjncy[ei].revPtr)->v;
    if(vi<vj) {
      adjncy[ei].mtd = 0;
      ntriangles = intersectLists(vault, modifiedGraph, hmap, vi, vj, (void*)&adjncy[ei], UPDATE_MTD);
      adjncy[ei].revPtr->mtd = adjncy[ei].mtd;
    }
  }

  for(int32_t i=0; i<HMAP_SIZE; i++) {
    edgeblock_t* ptr = hmap[i];
    while(ptr!=NULL) {
      vi = ptr->edge.v;
      vj = ptr->edge.revPtr->v;
      if(vi<vj) {
        ptr->edge.mtd = 0;
        ntriangles = intersectLists(vault, modifiedGraph, hmap, vi, vj, (void *)&(ptr->edge), UPDATE_MTD);
        ptr->edge.revPtr->mtd = ptr->edge.mtd;
      }
      ptr=ptr->next;
    }
  }
}

static int64_t updateTrussNumbers(vault_t *vault, ggraph_t *modifiedGraph, edgeblock_t **hmap, int32_t vi, int32_t vj, edgeblock_t* inserted_edge, rootEdges_struct *rootEdges) {
  /* Run the algorithm for each K(.) from local_ktmax down to 2 */
  int8_t isKFound = 0;
  int32_t suffixsum = 0;
  int64_t ntriangles = 0;
  int32_t ntriangles_update = 0;
  for(int i=vault->ktmax; i>rootEdges->ktmax; i--)
    suffixsum = suffixsum + rootEdges->queue[i].rear + 1;

  if(rootEdges->ktmax == -1) {
    inserted_edge->edge.k = 0;
    inserted_edge->edge.revPtr->k = 0;
    return 0;
  }

  int32_t visitedEdgeCount = 0;
  int32_t evictedEdgeCount = 0;
  for(int32_t k=rootEdges->ktmax; k>=0; k--) {
    edge_q *q = &(rootEdges->queue[k]);
    edge_q *searchQueue = gk_malloc(sizeof(edge_q), "updateTrussNumbers : Initial search space");
    searchSpace_struct search_s;
    search_s.queue = searchQueue;
    search_s.k = k;
    searchQueue->front = searchQueue->rear = -1;
    searchQueue->size = SEARCHSPACE_QUEUE_SIZE;
    searchQueue->edgelist = gk_malloc(sizeof(eedge_t*) * SEARCHSPACE_QUEUE_SIZE, "updateTrussNumbers : queue of edges in the initial search space");

    /* First insert the rootedges into the search space */
    int32_t f=q->front, r=q->rear;
    while(f != r) {
      eedge_t* cur_edge = q->edgelist[++f];
      cur_edge->visited = 1;
      cur_edge->revPtr->visited = 1;
      searchQueue->rear++;
      if(searchQueue->rear == searchQueue->size) {
        searchQueue->size *= 2;
        searchQueue->edgelist = gk_realloc(searchQueue->edgelist,  sizeof(eedge_t*) * searchQueue->size, "updateTrussNumbers : realloc queue of edges in the initial search space");
      }
      searchQueue->edgelist[searchQueue->rear] = cur_edge;
    }

    /*
    Explore the search space, while calculating the truss-degree.
    Insert the edges with insuffienct truss-degree into the eviction queue
    */
    edge_q *evictionList = gk_malloc(sizeof(edge_q), "updateTrussNumbers : List of edges int he search space that need to be evicted");
    evictionList->front = evictionList->rear = -1;
    evictionList->size = SEARCHSPACE_QUEUE_SIZE;
    evictionList->edgelist = gk_malloc(sizeof(eedge_t*) * SEARCHSPACE_QUEUE_SIZE, "updateTrussNumbers : queue of edges in the initial search space");
    f=searchQueue->front;
    while(f != searchQueue->rear) {
      eedge_t* cur_edge = searchQueue->edgelist[++f];
      int32_t Vj = cur_edge->v, Vi = cur_edge->revPtr->v;
      search_s.e = cur_edge;
      ntriangles += intersectLists(vault, modifiedGraph, hmap, Vi, Vj, (void *)&search_s, EXPLORE_SEARCHSPACE);
      if(cur_edge->td <= k) {
        evictionList->rear++;
        if(evictionList->rear == evictionList->size) {
          evictionList->size *= 2;
          evictionList->edgelist = gk_realloc(evictionList->edgelist,  sizeof(eedge_t*) * evictionList->size, "updateTrussNumbers : realloc queue of edges in the initial search space");
        }
        evictionList->edgelist[evictionList->rear] = cur_edge;
        cur_edge->visited = -1; /* visited = -1 corresponds to an edge being added to the eviction queue */
        cur_edge->revPtr->visited = -1;
      }
    }
    visitedEdgeCount += (searchQueue->rear+1);
    evictedEdgeCount += (evictionList->rear+1);

    /* Evict the edges in the eviction queue */
    search_s.queue = evictionList;
    f = evictionList->front;
    while(f != evictionList->rear) {
      eedge_t* cur_edge = evictionList->edgelist[++f];
      int32_t Vj = cur_edge->v, Vi = cur_edge->revPtr->v;
      search_s.e = cur_edge;
      ntriangles += intersectLists(vault, modifiedGraph, hmap, Vi, Vj, (void *)&search_s, EVICT_EDGES);
      cur_edge->visited = -2; /* visited = -2 corresponds to this edge being evicted,
                                 and all support triangles having their truss-degree reduced */
      cur_edge->revPtr->visited = -2;
    }

    /*
    If the truss number of the inserted edge is still not determined, then check if current 'k' satisfies
    Otherwise, just check if other edges that aren't evicted can have their K(.) increased
    */
    if(!isKFound) {
      /*
      K(.) of the inserted edge is not yet found. So, check if the count of the rootEdges
      that are not evicted provide enough support to the new edge to be a part of (k+1)-truss
      */
      int32_t f=q->front, cnt=0;
      while(f != q->rear) {
        eedge_t* cur_edge = q->edgelist[++f];
        if(cur_edge->visited==1)
          cnt++;
      }

      if(suffixsum+cnt > k) {
        isKFound = 1;
        inserted_edge->edge.k = k+1;
        inserted_edge->edge.revPtr->k = k+1;
        // inserted_edge->edge.mtd = suffixsum+cnt;
        // inserted_edge->edge.revPtr->mtd = suffixsum+cnt;
        // ntriangles_update += intersectLists(vault, modifiedGraph, hmap, inserted_edge->edge.v, inserted_edge->edge.revPtr->v, (void *)&(inserted_edge->edge.k), UPDATE_MTD);
        vault->ktmax = gk_max(vault->ktmax, k+1);
      }
      else {
        suffixsum = suffixsum + rootEdges->queue[k].rear + 1;
        f=searchQueue->front;
        while(f != searchQueue->rear) {
          eedge_t* cur_edge = searchQueue->edgelist[++f];
          cur_edge->visited=0;
          cur_edge->revPtr->visited=0;
          cur_edge->td=0;
          cur_edge->revPtr->td=0;
        }
      }
    }
    if(isKFound) {
      /* Increase the truss-numbers of the edges that were not evicted. */
      f=searchQueue->front;
      while(f != searchQueue->rear) {
        eedge_t* cur_edge = searchQueue->edgelist[++f];
        if(cur_edge->visited==1) {
          cur_edge->k+=1;
          cur_edge->revPtr->k+=1;
          // cur_edge->mtd = cur_edge->revPtr->mtd = cur_edge->td;
          // ntriangles_update += intersectLists(vault, modifiedGraph, hmap, cur_edge->v, cur_edge->revPtr->v, (void *)&cur_edge, UPDATE_MTD);
        }
        cur_edge->visited=0;
        cur_edge->revPtr->visited=0;
        cur_edge->td=0;
        cur_edge->revPtr->td=0;
      }
    }

    /* Free Resources */
    gk_free((void **)&searchQueue->edgelist, LTERM);
    gk_free((void **)&searchQueue, LTERM);
    gk_free((void **)&evictionList->edgelist, LTERM);
    gk_free((void **)&evictionList, LTERM);
  }
  printf("---------- Summary:-\nVisited Edges : %d; Evicted Edges : %d; Useful Work Done: %f \n----------\n", visitedEdgeCount, evictedEdgeCount, (float)(visitedEdgeCount-evictedEdgeCount)/visitedEdgeCount);
  printf("Additional triangles visited to update : %d\n", ntriangles_update);
  return ntriangles + ntriangles_update;
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

void stream_v4(params_t *params, vault_t *vault, char* edgesFile, char* outputLocation) {
  double totalRuntime = 0.0, addEdgeTime = 0.0, selectRootEdgesTime = 0.0, traversalTime = 0.0;
  int64_t ntriangles = 0;

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

    printf("\n**********************************\nNew Edge %d : (%d, %d)\n", idx+1, vault->iperm[vi]+1, vault->iperm[vj]+1);
    gk_startwctimer(vault->timer_1);
    /* 1. Add a new edge to the hashmap */
    edgeblock_t* inserted_edge = addNewEdge(vault, modifiedGraph, hmap, vi, vj);
    // printf("Added new edge\n");
    // printAdjList(modifiedGraph, hmap, vj);
    // printAdjList(modifiedGraph, hmap, vi);
    gk_stopwctimer(vault->timer_1);

    /* 2. Find root edges */
    gk_startwctimer(vault->timer_2);
    // printf("\n---> Select Root Edges\n");
    rootEdges_struct *rootEdges = selectRootEdges(vault, modifiedGraph, hmap, vi, vj);
    gk_stopwctimer(vault->timer_2);

    gk_startwctimer(vault->timer_3);
    // printf("\n---> Update\n");
    ntriangles += updateTrussNumbers(vault, modifiedGraph, hmap, vi, vj, inserted_edge, rootEdges);
    gk_stopwctimer(vault->timer_3);

    updateMTD(vault, modifiedGraph, hmap, inserted_edge);
    // ntriangles += intersectLists(vault, modifiedGraph, hmap, vi, vj, NULL, UPDATE_ROOT_EDGES);

    gk_stopwctimer(vault->timer_global);

    /* Write to output file */
#ifdef TEST_CORRECTNESS
    char* outputFile = gk_malloc(strlen(outputLocation) + 100, "stream");
    char* outputNum = gk_malloc(100, "stream");
    sprintf(outputNum, "%d", idx+1);
    strcpy(outputFile, outputLocation);
    strcat(outputFile, outputNum); strcat(outputFile, ".out");
    printf("Writing start...\n");
    writeKTToFile(vault, modifiedGraph, hmap, outputFile);
    printf("Writing done...\n");
    gk_free((void **)&outputFile, LTERM);
    gk_free((void **)&outputNum, LTERM);
#endif

    /* Output the runtimes to file */
    addEdgeTime = addEdgeTime + gk_getwctimer(vault->timer_1);
    selectRootEdgesTime = selectRootEdgesTime + gk_getwctimer(vault->timer_2);
    traversalTime = traversalTime + gk_getwctimer(vault->timer_3);
    totalRuntime = totalRuntime + gk_getwctimer(vault->timer_global);
    printf("Till now addEdgeTime Runtime : %.8lf\n", gk_getwctimer(vault->timer_1));
    printf("Till now selectRootEdgesTime Runtime : %.8lf\n", gk_getwctimer(vault->timer_2));
    printf("Till now Incremental traversalTime : %.8lf\n", gk_getwctimer(vault->timer_3));
    printf("Till now Incremental Runtime : %.8lf\n", gk_getwctimer(vault->timer_global));

    char* incrementalTimingsFile = params->infile;
    printf("Incremental Timings writtne to the file : %s \n", incrementalTimingsFile);
    char *ptr = strrchr(incrementalTimingsFile, '/');
    strcpy(ptr, "/incrementalTimings.txt");
    FILE *fpoutTimes;
    fpoutTimes = gk_fopen(incrementalTimingsFile, "a", "fpoutTimes");
    fprintf(fpoutTimes, "%.6lf\n", gk_getwctimer(vault->timer_global));
    gk_fclose(fpoutTimes);

    /* Free resources */
    for(int8_t i=0; i<=rootEdges->ktmax; i++)
      gk_free((void **)&rootEdges->queue[i].edgelist, LTERM);
    gk_free((void **)&rootEdges->queue, LTERM);
    gk_free((void **)&rootEdges, LTERM);
  }

  printf("\n\n\nTotal addEdgeTime Runtime : %.8lf\n", addEdgeTime);
  printf("Total selectRootEdgesTime Runtime : %.8lf\n", selectRootEdgesTime);
  printf("Total Incremental traversalTime : %.8lf\n", traversalTime);
  printf("Total Incremental Runtime : %.8lf\n", totalRuntime);
  printf("Total Incremental Ops : %"PRId64"\n", ntriangles);

  /* Cleanup */
  gk_free((void **)&newEdges, LTERM);
  gk_free((void **)&vault->kt, LTERM);
  gk_free((void **)&vault->endIdx, LTERM);
  gk_free((void **)&vault->revIdx, LTERM);
}
