/*!
\file
\brief The various k-truss decomposition routines
\date Started 6/3/2017
\author George
\version\verbatim $Id: cmdline.c 20946 2017-05-10 23:12:48Z karypis $ \endverbatim
*/

#include "kt.h"

#define hfun1(vi, vj, i, range) \
           (((ssize_t)(((((ssize_t)vi)+5)^((ssize_t)vj)*(((ssize_t)vi>>32)+1)^((ssize_t)vj<<7)) + (i)*((1+((ssize_t)vi>>3)+1)^((ssize_t)vj<<5))))%range)


#ifndef DYNAMIC_CHUNK
#define DYNAMIC_CHUNK 16
#endif

/*************************************************************************/
/*! Determine the iperm for the key order using counting sort.
*/
/*************************************************************************/
int32_t *gk_i32kvipermi(int32_t n, gk_i32kv_t *cand)
{
  int i, j, k, range;
  int32_t *counts, *iperm;

  for (range=0, i=0; i<n; i++) {
    if (cand[i].key > range)
      range = cand[i].key;
  }
  range++;

  counts = gk_i32smalloc(range+1, 0, "counts");
  for (i=0; i<n; i++)
    counts[cand[i].key]++;
  MAKECSR(i, range, counts);

  iperm = gk_i32smalloc(n, 0, "iperm");
  for (i=0; i<n; i++)
    iperm[counts[cand[i].key]++] = i;

  gk_free((void **)&counts, LTERM);

  return iperm;
}


/*************************************************************************/
/*! Reorder the vertices in the graph in inc degree and return the upper
    triangular part of the reordered graph in which the adjancency lists
    are sorted in increasing order.
*/
/*************************************************************************/
gk_graph_t *kt_PreprocessAndExtractUpper(params_t *params, vault_t *vault)
{
  int32_t vi, vj, vk, nvtxs;
  ssize_t ei, eiend, ej, ejend, nedges;
  ssize_t *xadj, *uxadj;
  int32_t *adjncy, *uadjncy, *perm=NULL, *iperm=NULL;
  gk_i32kv_t *cand=NULL;
  gk_graph_t *graph;

  nvtxs  = vault->graph->nvtxs;
  xadj   = vault->graph->xadj;
  adjncy = vault->graph->adjncy;

  cand = gk_i32kvmalloc(nvtxs, "cand");
  for (vi=0; vi<nvtxs; vi++) {
    cand[vi].key = (int32_t)(xadj[vi+1]-xadj[vi]);
    cand[vi].val = vi;
  }

  perm  = vault->perm  = gk_i32smalloc(nvtxs, -1, "perm");   /* perm[old-vtx-num]  => new-vtx-num */
  iperm = vault->iperm = gk_i32kvipermi(nvtxs, cand);        /* iperm[new-vtx-num] => old-vtx-num */
  for (vi=0; vi<nvtxs; vi++)
    perm[iperm[vi]] = vi;


  /* create the reordered/sorted upper triangular portion of the graph */
  graph = gk_graph_Create();
  graph->nvtxs  = nvtxs;
  graph->xadj   = uxadj   = gk_zmalloc(nvtxs+1, "uxadj");
  graph->adjncy = uadjncy = gk_i32malloc(10+(xadj[nvtxs]>>1), "uadjncy");

  uxadj[0] = nedges = 0;
  for (vi=0; vi<nvtxs; vi++) {
    vj = iperm[vi];
    for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
      assert(adjncy[ej] < nvtxs);
      if ((vk = perm[adjncy[ej]]) > vi) /* keep only the upper part */
        uadjncy[nedges++] = vk;
    }
    uxadj[vi+1] = nedges;

    if (nedges-uxadj[vi] > 1)
      gk_i32sorti(nedges-uxadj[vi], uadjncy+uxadj[vi]);  /* sort adjncy list */
  }
  printf("Upper nedges: %zd out of %zd\n", uxadj[nvtxs], xadj[nvtxs]);

  gk_free((void **)&cand, LTERM);

  return graph;
}


/*************************************************************************/
/*! Creates the transpose of the upper-triangular graph with location
    offsets at +1 locations.
    This is used for the JIK algorithm.
*/
/*************************************************************************/
gk_graph_t *kt_TransposeUforJIK(params_t *params, gk_graph_t *graph)
{
  int32_t vi, vj, nvtxs;
  ssize_t ei, eiend, nedges;
  ssize_t *xadj, *txadj;
  int32_t *adjncy, *tadjncy;
  gk_graph_t *tgraph;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  tgraph         = gk_graph_Create();
  tgraph->nvtxs  = nvtxs;
  tgraph->xadj   = txadj   = gk_zsmalloc(nvtxs+1, 0, "txadj");
  tgraph->adjncy = tadjncy = gk_i32malloc(2*(xadj[nvtxs]+1), "tadjncy");

  for (vi=0; vi<nvtxs; vi++) {
    if (xadj[vi+1]-xadj[vi] < 2)
      continue;
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend-1; ei++)
      txadj[adjncy[ei]] += 2;
  }
  MAKECSR(vi, nvtxs, txadj);

  for (vi=0; vi<nvtxs; vi++) {
    if (xadj[vi+1]-xadj[vi] < 2)
      continue;

    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend-1; ei++) {
      vj = adjncy[ei];
      tadjncy[txadj[vj]++] = vi;
      tadjncy[txadj[vj]++] = ei-xadj[vi]+1;  /* row-offset */
    }
  }
  SHIFTCSR(vi, nvtxs, txadj);

  return tgraph;
}


/*************************************************************************/
/*! Checks if the supports computed by the TC code is correct.
*/
/*************************************************************************/
void kt_CheckInitialSupport(params_t *params, vault_t *vault)
{
  int32_t uvi, vi, vik, vj, vjk, vk, nvtxs, nh;
  ssize_t uei, ei, ej;
  ssize_t *xadj, *uxadj;
  int32_t *adjncy, *uadjncy, *uadjwgt;
  int32_t *map;

  nvtxs  = vault->graph->nvtxs;
  xadj   = vault->graph->xadj;
  adjncy = vault->graph->adjncy;

  uxadj   = vault->ugraph->xadj;
  uadjncy = vault->ugraph->adjncy;
  uadjwgt = vault->ugraph->iadjwgt;

  map = gk_i32smalloc(nvtxs, -1, "map");

  for (uvi=0; uvi<nvtxs; uvi++) {
    vi = vault->iperm[uvi];
    for (ei=xadj[vi]; ei<xadj[vi+1]; ei++)
      map[adjncy[ei]] = vi;

    for (uei=uxadj[uvi]; uei<uxadj[uvi+1]; uei++) {
      vj = vault->iperm[uadjncy[uei]];
      nh = uadjwgt[uei];
      for (ej=xadj[vj]; ej<xadj[vj+1]; ej++)
        nh -= (map[adjncy[ej]] == vi ? 1 : 0);

      GKASSERT(nh == 0);
    }
  }

  gk_free((void **)&map, LTERM);
}

/*************************************************************************/
/*! Checks if the supports computed by the TC code is correct.
*/
/*************************************************************************/
void kt_CheckKTrussDecomposition(params_t *params, vault_t *vault)
{
  int32_t k, vi, vj, vk, nvtxs, knvtxs, nh;
  ssize_t ei, ej, nedges;
  ssize_t *xadj;
  int32_t *adjncy;
  int32_t *map;
  ktedge_t *ktedges;

  nvtxs   = vault->graph->nvtxs;
  nedges  = vault->nedges;
  ktedges = vault->ktedges;

  for (k=1; k<=vault->ktmax; k++) {
    xadj = gk_zsmalloc(nvtxs+1, 0, "xadj");
    for (ei=0; ei<nedges; ei++) {
      if (ktedges[ei].k >= k) {
        xadj[ktedges[ei].vi]++;
        xadj[ktedges[ei].vj]++;
      }
    }
    for (knvtxs=0, vi=0; vi<nvtxs; vi++)
      knvtxs += (xadj[vi] > 0 ? 1 : 0);
    MAKECSR(vi, nvtxs, xadj);

    adjncy = gk_i32malloc(xadj[nvtxs], "adjncy");
    for (ei=0; ei<nedges; ei++) {
      if (ktedges[ei].k >= k) {
        adjncy[xadj[ktedges[ei].vi]++] = ktedges[ei].vj;
        adjncy[xadj[ktedges[ei].vj]++] = ktedges[ei].vi;
      }
    }
    SHIFTCSR(vi, nvtxs, xadj);

    map = gk_i32smalloc(nvtxs, -1, "map");
    for (vi=0; vi<nvtxs; vi++) {
      for (ei=xadj[vi]; ei<xadj[vi+1]; ei++)
        map[adjncy[ei]] = vi;

      for (ei=xadj[vi]; ei<xadj[vi+1]; ei++) {
        vj = adjncy[ei];
        for (nh=0, ej=xadj[vj]; ej<xadj[vj+1]; ej++)
          nh += (map[adjncy[ej]] == vi ? 1 : 0);

        GKASSERT(nh >= k);
      }
    }

    printf("k-truss: %4d, nvtxs: %7d, nedges: %8zd\n", k+2, knvtxs, xadj[nvtxs]);

    gk_free((void **)&xadj, &adjncy, &map, LTERM);
  }
}


/*************************************************************************/
/*! Takes the sups[] array associated with the edges and creates the
    ktedges information in the vault.
*/
/*************************************************************************/
void kt_Sups2KTEdges(params_t *params, vault_t *vault, int32_t ktmax, int32_t *sups)
{
  int32_t vi, nvtxs;
  ssize_t ei, eiend, ej, nedges;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;

  if (params->outfile == NULL)
    return;

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;
  adjwgt = vault->ugraph->iadjwgt;

  vault->nedges = xadj[nvtxs];
  vault->ktmax  = ktmax;
  vault->ktedges = (ktedge_t *)gk_malloc(xadj[nvtxs]*sizeof(ktedge_t), "ktedges");

  for (ej=0, nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++, ej++) {
      vault->ktedges[ej].vi = gk_min(vault->iperm[vi], vault->iperm[adjncy[ei]]);
      vault->ktedges[ej].vj = gk_max(vault->iperm[vi], vault->iperm[adjncy[ei]]);
      if (adjwgt[ei] > 0)
        vault->ktedges[ej].k = -sups[nedges++] + 2;
      else
        vault->ktedges[ej].k = 2;
    }
  }
}

/*************************************************************************/
/*! The hash-map-based edge-triangle-support counting routine that uses
    the JIK triangle enumeration scheme.

    This is the mapjikv2 tc version.
*/
/*************************************************************************/
int64_t kt_ComputeEdgeSupport(params_t *params, vault_t *vault)
{
  int32_t vi, vj, vk, vl, nvtxs, nlocal;
  ssize_t ei, eiend, ej, ejstart, ejend;
  int64_t ntriangles, ntriangles2;
  ssize_t *xadj, *txadj;
  int32_t *adjncy, *tadjncy, *adjwgt;
  int32_t l, tnc, nc, hmsize, tlsize, tlstart, *hmap, *tmap;

  gk_startwctimer(vault->timer_2);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;
  adjwgt = vault->ugraph->iadjwgt;

  txadj   = vault->lgraph->xadj;
  tadjncy = vault->lgraph->adjncy;


  /* determine the size of the hash-map and convert it into a format
     that is compatible with a bitwise AND operation */
  for (hmsize=0, vi=0; vi<nvtxs; vi++)
    hmsize = gk_max(hmsize, (int32_t)(xadj[vi+1]-xadj[vi]));
  for (l=1; hmsize>(1<<l); l++);
  hmsize = (1<<(l+4))-1;
  hmap = gk_i32smalloc(hmsize+1, -1, "hmap");
  printf("& compatible maximum hmsize: %"PRId32"\n", hmsize);

  /* determine the size of the tail-map and allocate memory for it */
  for (vi=(nvtxs>>2); vi<nvtxs; vi++) {
    if ((txadj[vi+1]-txadj[vi])<<9 > vi)
      break;
    if ((xadj[vi+1]-xadj[vi])<<4 > nvtxs-vi)
      break;
  }
  tlsize  = nvtxs - vi + 100;
  tlstart = nvtxs-tlsize;
  tmap    = gk_i32smalloc(tlsize, -1, "tmap");
  tmap   -= tlstart; /* make indexing simpler */
  printf("tlsize: %"PRId32"\n", tlsize);


  /* start counting triangles */
  if (params->dbglvl&1)
    gk_startwctimer(vault->timer_4);

  /* use a combination of hmap and tmap */
  ntriangles = 0;
  hmsize     = 0;
  tnc        = 0;
  for (vj=1; vj<tlstart; vj++) {
    if (xadj[vj+1] == xadj[vj] || txadj[vj+1] == txadj[vj])
      continue;

    /* if needed, increase the working hmsize */
    if ((xadj[vj+1]-xadj[vj])<<3 > 1 + (hmsize>>4) + (hmsize>>1)) {
      hmsize = xadj[vj+1]-xadj[vj];
      for (l=1; hmsize>(1<<l); l++);
      hmsize = (1<<(l+4))-1;

      if (params->dbglvl&1) {
        gk_stopwctimer(vault->timer_4);
        printf("vj: %9d tlstart: %d degree: %5zd %7zd hmsize: %6d tnc: %7d time: %5.2lfs\n",
            vj, tlstart, xadj[vj+1]-xadj[vj], txadj[vj+1]-txadj[vj],
            hmsize, tnc, gk_getwctimer(vault->timer_4));
        tnc = 0;
        gk_clearwctimer(vault->timer_4);
        gk_startwctimer(vault->timer_4);
      }
    }

    /* hash Adj(vj) using hmap for the front and tmap for the last tlsize indices */
    for (nc=0, ej=ejstart=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
      if ((vk = adjncy[ej]) >= tlstart)
        break;
      for (l=(vk&hmsize); hmap[l]!=-1; l=((l+1)&hmsize), nc++);
      hmap[l] = ej-ejstart;
    }
    for (; ej<ejend; ej++)
      tmap[adjncy[ej]] = ej-ejstart;

    tnc += nc;

    /* find intersections */
    if (nc > 0) { /* we had collisions */
      for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
        vi = tadjncy[ej];
        for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
          if ((vk = adjncy[ei]) >= tlstart)
            break;
          l = vk&hmsize;
          if (hmap[l] == -1)
            continue;
          if (adjncy[ejstart+hmap[l]] == vk) {
            adjwgt[ei]++;
            adjwgt[ejstart+hmap[l]]++;
            nlocal++;
            continue;
          }
          for (l=((l+1)&hmsize); hmap[l]!=-1 && adjncy[ejstart+hmap[l]]!=vk; l=((l+1)&hmsize));
          if (hmap[l]!=-1 && adjncy[ejstart+hmap[l]] == vk) {
            adjwgt[ei]++;
            adjwgt[ejstart+hmap[l]]++;
            nlocal++;
          }
        }
        for (; ei<eiend; ei++) {
          if (tmap[adjncy[ei]] != -1) {
            assert(adjncy[ejstart+tmap[adjncy[ei]]] == adjncy[ei]);
            adjwgt[ei]++;
            adjwgt[ejstart+tmap[adjncy[ei]]]++;
            nlocal++;
          }
        }

        if (nlocal > 0) {
          ntriangles += nlocal;

          assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
          adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
        }
      }

      /* reset hmap/tmap */
      for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
        if ((vk = adjncy[ej]) >= tlstart)
          break;
        for (l=(vk&hmsize); hmap[l]==-1 || adjncy[ejstart+hmap[l]]!=vk; l=((l+1)&hmsize));
        hmap[l] = -1;
      }
      for (; ej<ejend; ej++)
        tmap[adjncy[ej]] = -1;
    }
    else { /* there were no collisons */
      for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
        vi = tadjncy[ej];

        for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
          if ((vk = adjncy[ei]) >= tlstart)
            break;
          if (hmap[vk&hmsize]!=-1 && adjncy[ejstart+hmap[vk&hmsize]] == vk) {
            adjwgt[ei]++;
            adjwgt[ejstart+hmap[vk&hmsize]]++;
            nlocal++;
          }
        }
        for (; ei<eiend; ei++) {
          if (tmap[adjncy[ei]] != -1) {
            assert(adjncy[ejstart+tmap[adjncy[ei]]] == adjncy[ei]);
            adjwgt[ei]++;
            adjwgt[ejstart+tmap[adjncy[ei]]]++;
            nlocal++;
          }
        }

        if (nlocal > 0) {
          ntriangles += nlocal;

          assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
          adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
        }
      }

      /* reset hmap/tmap */
      for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
        if ((vk = adjncy[ej]) >= tlstart)
          break;
        hmap[vk&hmsize] = -1;
      }
      for (; ej<ejend; ej++)
        tmap[adjncy[ej]] = -1;
    }
  }

  printf("ntriangles: %zd\n", ntriangles);

  if (params->dbglvl&1) {
    gk_stopwctimer(vault->timer_4);
    printf("vj: %9d tlstart: %d degree: %5zd %7zd hmsize: %6d tnc: %7d time: %5.2lfs\n",
        vj, tlstart, xadj[vj+1]-xadj[vj], txadj[vj+1]-txadj[vj],
        hmsize, tnc, gk_getwctimer(vault->timer_4));
    tnc = 0;
    gk_clearwctimer(vault->timer_4);
    gk_startwctimer(vault->timer_4);
  }

  /* use tmap for the last tlsize rows */
  for (; vj<nvtxs; vj++) {
    if (1 || xadj[vj+1]-xadj[vj] < nvtxs-vj-1) {
      /* hash Adj(vj) */
      for (ej=ejstart=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++)
        tmap[adjncy[ej]] = ej-ejstart;

      /* find intersections */
      for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
        vi = tadjncy[ej];
        for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
          if (tmap[adjncy[ei]] != -1) {
            adjwgt[ei]++;
            adjwgt[ejstart+tmap[adjncy[ei]]]++;
            nlocal++;
          }
        }

        if (nlocal > 0) {
          ntriangles += nlocal;

          assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
          adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
        }
      }

      /* reset tmap */
      for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++)
        tmap[adjncy[ej]] = -1;
    }
    else { /* the row is dense */  /* TODO: This has not been updated */
      tnc++;
      /* find intersections */
      for (nlocal=0, ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
        vi = tadjncy[ej];
        nlocal += xadj[vi+1]-xadj[vi]-tadjncy[ej+1];
      }
      ntriangles += nlocal;
    }
  }
  gk_stopwctimer(vault->timer_2);

  if (params->dbglvl&1) {
    gk_stopwctimer(vault->timer_4);

    vj = nvtxs-2;
    printf("vj: %9d tlstart: %d degree: %5zd %7zd hmsize: %6d tnc: %7d time: %5.2lfs\n",
        vj, tlstart, xadj[vj+1]-xadj[vj], txadj[vj+1]-txadj[vj],
        hmsize, tnc, gk_getwctimer(vault->timer_4));

    for (ntriangles2=0, ei=0; ei<xadj[nvtxs]; ei++)
      ntriangles2 += adjwgt[ei];

    printf("Sanity check: ntriangles: %zd %zd %zd\n", ntriangles, ntriangles2/3, ntriangles2%3);
  }

  tmap += tlstart;
  gk_free((void **)&hmap, &tmap, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline1(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj, sup;
  } *edges;

  struct aii_s {
    int32_t vj, inc;
    size_t id;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;


  gk_ipq_t *queue;

  int32_t k, vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eip, eiend, ej, ejp, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the increment (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi  = vi;
        edges[nedges].vj  = adjncy[ei];
        edges[nedges].sup = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].id  = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].id  = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  gk_startwctimer(vault->timer_1);
  /* create and populate the priority queue
     [the negative on the key is because my priority queues are max priority queues */
  queue = gk_ipqCreate(nedges);
  for (ei=0; ei<nedges; ei++)
    gk_ipqInsert(queue, ei, -edges[ei].sup);
  gk_stopwctimer(vault->timer_1);

  printf("Pqueue stats: nnodes: %zu, nedges: %zd\n", gk_ipqLength(queue), nedges);
  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1;;k++) {
    if (gk_ipqLength(queue) == 0)
      break;

    nlocal = 0;
    nltriangles = 0;
    while (k >= -gk_ipqSeeTopKey(queue)) {
      if (gk_ipqLength(queue) == 0)
        break;

      nlocal++;

      ti = gk_ipqGetTop(queue);
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (edges[ti].sup > 0) {
        gk_startwctimer(vault->timer_3);
        sup = edges[ti].sup;
        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, edges[ti].sup,
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        ei    = xaii[vi].start;
        eiend = xaii[vi].end;
        vik   = aii[ei].vj;

        ej    = xaii[vj].start;
        ejend = xaii[vj].end;
        vjk   = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei < eiend && ej < ejend) {
          //printf("  vik: %d  vjk: %d\n", vik, vjk);
          if (vik < vjk) {
            ei  += aii[ei].inc;
            vik  = aii[ei].vj;
          }
          else if (vjk < vik) {
            ej  += aii[ej].inc;
            vjk  = aii[ej].vj;
          }
          else {
            //printf("  ** vk: %d\n", vik);
            nltriangles++;

            edges[aii[ei].id].sup--;
            gk_ipqUpdate(queue, aii[ei].id, -edges[aii[ei].id].sup);
            ei  += aii[ei].inc;
            vik  = aii[ei].vj;

            edges[aii[ej].id].sup--;
            gk_ipqUpdate(queue, aii[ej].id, -edges[aii[ej].id].sup);
            ej  += aii[ej].inc;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        gk_stopwctimer(vault->timer_3);
      }

      /* remove the edge from both adjacency lists */
      for (eip=-1, ei=xaii[vi].start; aii[ei].vj<vj; eip=ei, ei+=aii[ei].inc);
      ASSERT(aii[ei].vj == vj);
      if (eip == -1)
        xaii[vi].start += aii[ei].inc;
      else
        aii[eip].inc += aii[ei].inc;

      edges[aii[ei].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */

      for (ejp=-1, ej=xaii[vj].start; aii[ej].vj<vi; ejp=ej, ej+=aii[ej].inc);
      ASSERT(aii[ej].vj == vi);
      if (ejp == -1)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ejp].inc += aii[ej].inc;

      edges[aii[ej].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (nlocal > 0 && params->dbglvl&1)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  gk_free((void **)&edges, &aii, &xaii, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline2(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj, sup;
  } *edges;

  struct aii_s {
    int32_t vj, inc;
    size_t id;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  gk_ipq_t *queue;

  int32_t k, vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eip, eiend, ej, ejp, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  ssize_t nupdates, *updindices, *updmark;

  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the increment (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi  = vi;
        edges[nedges].vj  = adjncy[ei];
        edges[nedges].sup = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].id  = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].id  = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  gk_startwctimer(vault->timer_1);
  /* create and populate the priority queue
     [the negative on the key is because my priority queues are max priority queues */
  queue = gk_ipqCreate(nedges);
  for (ei=0; ei<nedges; ei++)
    gk_ipqInsert(queue, ei, -edges[ei].sup);
  gk_stopwctimer(vault->timer_1);

  printf("Pqueue stats: nnodes: %zu, nedges: %zd\n", gk_ipqLength(queue), nedges);

  updindices = gk_zmalloc(nedges, "updindices");
  updmark    = gk_zsmalloc(nedges, -1, "updmark");
  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1;;k++) {
    if (gk_ipqLength(queue) == 0)
      break;

    nlocal = 0;
    nltriangles = 0;
BACK:
    nupdates = 0;
    while (k >= -gk_ipqSeeTopKey(queue)) {
      if (gk_ipqLength(queue) == 0)
        break;

      nlocal++;

      ti = gk_ipqGetTop(queue);
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (edges[ti].sup > 0) {
        gk_startwctimer(vault->timer_3);
        sup = edges[ti].sup;
        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, edges[ti].sup,
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        ei    = xaii[vi].start;
        eiend = xaii[vi].end;
        vik   = aii[ei].vj;

        ej    = xaii[vj].start;
        ejend = xaii[vj].end;
        vjk   = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei < eiend && ej < ejend) {
          //printf("  vik: %d  vjk: %d\n", vik, vjk);
          if (vik < vjk) {
            ei  += aii[ei].inc;
            vik  = aii[ei].vj;
          }
          else if (vjk < vik) {
            ej  += aii[ej].inc;
            vjk  = aii[ej].vj;
          }
          else {
            //printf("  ** vk: %d\n", vik);
            nltriangles++;

            edges[aii[ei].id].sup--;
            if (updmark[aii[ei].id] == -1 ) {
              updindices[nupdates++] = aii[ei].id;
              updmark[aii[ei].id] = 1;
            }
            ei  += aii[ei].inc;
            vik  = aii[ei].vj;

            edges[aii[ej].id].sup--;
            if (updmark[aii[ej].id] == -1 ) {
              updindices[nupdates++] = aii[ej].id;
              updmark[aii[ej].id] = 1;
            }
            ej  += aii[ej].inc;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        gk_stopwctimer(vault->timer_3);
      }

      /* remove the edge from both adjacency lists */
      gk_startwctimer(vault->timer_5);
      for (eip=-1, ei=xaii[vi].start; aii[ei].vj<vj; eip=ei, ei+=aii[ei].inc);
      ASSERT(aii[ei].vj == vj);
      if (eip == -1)
        xaii[vi].start += aii[ei].inc;
      else
        aii[eip].inc += aii[ei].inc;

      edges[aii[ei].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */
      gk_stopwctimer(vault->timer_5);

      gk_startwctimer(vault->timer_6);
      for (ejp=-1, ej=xaii[vj].start; aii[ej].vj<vi; ejp=ej, ej+=aii[ej].inc);
      ASSERT(aii[ej].vj == vi);
      if (ejp == -1)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ejp].inc += aii[ej].inc;

      edges[aii[ej].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */
      gk_stopwctimer(vault->timer_6);
    }

    if (nupdates > 0) {
      //printf("  nupdates: %zd [%zd]\n", nupdates, nltriangles);
      for (ei=0; ei<nupdates; ei++) {
        if (edges[updindices[ei]].sup >= 0)
          gk_ipqUpdate(queue, updindices[ei], -edges[updindices[ei]].sup);
        updmark[updindices[ei]] = -1;
      }
      goto BACK;
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  gk_free((void **)&edges, &aii, &xaii, &updindices, &updmark, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline3(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj, sup;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
    size_t id;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  gk_ipq_t *queue;

  int32_t k, vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eistart, eiend, ej, ejstart, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  ssize_t nupdates, *updindices, *updmark;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi  = vi;
        edges[nedges].vj  = adjncy[ei];
        edges[nedges].sup = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        aii[xaii[vj].start].id  = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        aii[xaii[vi].start].id  = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  gk_startwctimer(vault->timer_1);
  /* create and populate the priority queue
     [the negative on the key is because my priority queues are max priority queues */
  queue = gk_ipqCreate(nedges);
  for (ei=0; ei<nedges; ei++)
    gk_ipqInsert(queue, ei, -edges[ei].sup);
  gk_stopwctimer(vault->timer_1);

  printf("Pqueue stats: nnodes: %zu, nedges: %zd\n", gk_ipqLength(queue), nedges);

  updindices = gk_zmalloc(nedges, "updindices");
  updmark    = gk_zsmalloc(nedges, -1, "updmark");

  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1;;k++) {
    if (gk_ipqLength(queue) == 0)
      break;

    nltriangles = 0;
    nlocal = 0;
BACK:
    nupdates = 0;
    while (k >= -gk_ipqSeeTopKey(queue)) {
      if (gk_ipqLength(queue) == 0)
        break;

      nlocal++;

      ti = gk_ipqGetTop(queue);
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (edges[ti].sup > 0) {
        gk_startwctimer(vault->timer_3);
        sup = edges[ti].sup;
        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, edges[ti].sup,
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        ei      = xaii[vi].end-1;
        eistart = xaii[vi].start;
        vik     = aii[ei].vj;

        ej      = xaii[vj].end-1;
        ejstart = xaii[vj].start;
        vjk     = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei >= eistart && ej >= ejstart) {
          //printf("  vik: %d  vjk: %d\n", vik, vjk);
          if (vik > vjk) {
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;
          }
          else if (vjk > vik) {
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;
          }
          else {
            //printf("  ** vk: %d\n", vik);
            nltriangles++;

            edges[aii[ei].id].sup--;
            if (updmark[aii[ei].id] == -1 ) {
              updindices[nupdates++] = aii[ei].id;
              updmark[aii[ei].id] = 1;
            }
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;

            edges[aii[ej].id].sup--;
            if (updmark[aii[ej].id] == -1 ) {
              updindices[nupdates++] = aii[ej].id;
              updmark[aii[ej].id] = 1;
            }
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        gk_stopwctimer(vault->timer_3);
      }

      /* remove the edge from both adjacency lists */
      gk_startwctimer(vault->timer_5);
      for (ei=xaii[vi].start; aii[ei].vj<vj; ei+=aii[ei].inc);
      ASSERT(aii[ei].vj == vj);
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      edges[aii[ei].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */
      gk_stopwctimer(vault->timer_5);

      gk_startwctimer(vault->timer_6);
      for (ej=xaii[vj].start; aii[ej].vj<vi; ej+=aii[ej].inc);
      ASSERT(aii[ej].vj == vi);
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;

      edges[aii[ej].id].sup = -k;  /* this is used for encoding the maximal value of k of that edge */
      gk_stopwctimer(vault->timer_6);

    }

    if (nupdates > 0) {
      for (ei=0; ei<nupdates; ei++) {
        if (edges[updindices[ei]].sup >= 0)
          gk_ipqUpdate(queue, updindices[ei], -edges[updindices[ei]].sup);
        updmark[updindices[ei]] = -1;
      }
      goto BACK;
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  gk_free((void **)&edges, &aii, &xaii, &updindices, &updmark, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    list-traversal + priority queue + O(1) deletion
*/
/*************************************************************************/
int64_t kt_Baseline4(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  gk_ipq_t *queue;

  int32_t k, vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eistart, eiend, ej, ejstart, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  ssize_t nupdates, nmaxupdates, *updindices;
  int32_t *sups;
  ssize_t *ids;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);

  if (params->dbglvl&8)
    kt_CheckInitialSupport(params, vault);

  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii   = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii    = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges  = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups   = gk_i32malloc(nedges, "sups");
  ids    = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* create and populate the priority queue
     [the negative on the key is because my priority queues are max priority queues */
  queue = gk_ipqCreate(nedges);
  for (ei=0; ei<nedges; ei++)
    gk_ipqInsert(queue, ei, -sups[ei]);

  printf("Pqueue stats: nnodes: %zu, nedges: %zd\n", gk_ipqLength(queue), nedges);

  nmaxupdates = nedges + 2*nvtxs;
  updindices = gk_zmalloc(nmaxupdates, "updindices");

  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1;;k++) {
    if (gk_ipqLength(queue) == 0)
      break;

    nltriangles = 0;
    nlocal = 0;
BACK:
    nupdates = 0;
    while (k >= -gk_ipqSeeTopKey(queue)) {
      if (nupdates + 2*nvtxs > nmaxupdates)
        break;

      if (gk_ipqLength(queue) == 0)
        break;

      nlocal++;

      ti = gk_ipqGetTop(queue);
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      /* decrease the support of the intersection of the adjacency lists */
      if (sups[ti] > 0) {
        gk_startwctimer(vault->timer_3);
        sup = sups[ti];

        nltriangles += sup;

        ei      = xaii[vi].end-1;
        eistart = xaii[vi].start;
        vik     = aii[ei].vj;

        ej      = xaii[vj].end-1;
        ejstart = xaii[vj].start;
        vjk     = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei >= eistart && ej >= ejstart) {
          if (vik > vjk) {
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;
          }
          else if (vjk > vik) {
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;
          }
          else {
            updindices[nupdates++] = ids[ei];
            updindices[nupdates++] = ids[ej];

            sups[ids[ei]]--;
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;

            sups[ids[ej]]--;
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);
        gk_stopwctimer(vault->timer_3);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (nupdates > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdates; ei++) {
        if (sups[updindices[ei]] >= 0) {
          if (queue->heap[queue->locator[updindices[ei]]].key != -sups[updindices[ei]])
            gk_ipqUpdate(queue, updindices[ei], -sups[updindices[ei]]);
        }
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_free((void **)&edges, &aii, &xaii, &updindices, &sups, &ids, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    list-traversal + bucket-based supports + O(1) deletion
*/
/*************************************************************************/
int64_t kt_Baseline5(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nltriangles, sup;
  ssize_t ti, ei, eistart, eiend, ej, ejstart, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  int32_t k, sk, nsups, *sups;
  ssize_t si, eid, *ids, *sptr, *sind, *sloc;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  sptr = gk_zsmalloc(nsups+1, 0, "sptr"); /* the starts of the different support buckets */
  sind = gk_zmalloc(nedges, "sind");      /* the set of edge IDs in each support bucket */
  sloc = gk_zmalloc(nedges, "sloc");      /* the sloc[edge ID] => location in sind[] mapping */

  for (ei=0; ei<nedges; ei++)
    sptr[sups[ei]]++;
  MAKECSR(vi, nsups, sptr);
  for (ei=0; ei<nedges; ei++) {
    sloc[ei] = sptr[sups[ei]];
    sind[sptr[sups[ei]]++] = ei;
  }
  SHIFTCSR(vi, nsups, sptr);

  gk_stopwctimer(vault->timer_ktsetup);


  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && sptr[k]!=nedges; k++) {
    nltriangles = 0;
    for (si=sptr[k]; si<sptr[k+1]; si++) {
      ti = sind[si];
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        ei      = xaii[vi].end-1;
        eistart = xaii[vi].start;
        vik     = aii[ei].vj;

        ej      = xaii[vj].end-1;
        ejstart = xaii[vj].start;
        vjk     = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei >= eistart && ej >= ejstart) {
          if (vik > vjk) {
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;
          }
          else if (vjk > vik) {
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;
          }
          else {
            eid = ids[ei];
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;

            eid = ids[ej];
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (params->dbglvl&1)
      printf("k: %7d; nleft: %7zd, nltriangles: %7d\n", k, nedges-sptr[k+1], nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &sptr, &sind, &sloc, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    list-traversal + bucket-based supports + O(1) deletion
*/
/*************************************************************************/
int64_t kt_Baseline5b(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  struct slist_s {
    ssize_t neid, peid;
    int32_t sup;
  } *slist;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nltriangles, sup;
  ssize_t ti, ei, eistart, eiend, ej, ejstart, ejend;
  int64_t nedges, nleft, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  int32_t k, nsups, *sups;
  ssize_t *ids, *shead;
  ssize_t nupdates, nmaxupdates, *updindices;

  double timer_currk = 0.;

  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);

#if VERBOSE
  printf("supports:\n");
  for(int v=0; v < nvtxs; ++v) {
    for(ssize_t e=xadj[v]; e < xadj[v+1]; ++e) {
      printf("(%2d, %2d) perm[%2d, %2d] = %d\n",
          v+1, adjncy[e]+1,
          vault->iperm[v]+1, vault->iperm[adjncy[e]]+1,
          adjwgt[e]);
    }
  }
#endif

  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  /* the heads and "link list" that form the support buckets */
  shead = gk_zsmalloc(nsups, -1, "shead");
  slist = (struct slist_s *)gk_malloc((nedges+1)*sizeof(struct slist_s), "slist");

  slist++;  /* this is to allow slist[-1] to be valid */
  for (ei=0; ei<nedges; ei++) {
    slist[ei].sup  = sups[ei];
    slist[ei].peid = -1;
    slist[ei].neid = shead[sups[ei]];
    if (shead[sups[ei]] != -1)
      slist[shead[sups[ei]]].peid = ei;
    shead[sups[ei]] = ei;
  }

  nmaxupdates = nedges + 2*nvtxs;
  updindices = gk_zmalloc(nmaxupdates, "updindices");

  gk_stopwctimer(vault->timer_ktsetup);


  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;
  nleft      = nedges;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && nleft>0; k++) {
    nltriangles = 0;
    gk_clearwctimer(timer_currk);
    gk_startwctimer(timer_currk);

BACK:
    nupdates = 0;
    for (ti=shead[k]; ti!=-1; ti=slist[ti].neid) {
      if (nupdates + 2*nvtxs > nmaxupdates)
        break;

      nleft--;

      vi = edges[ti].vi;
      vj = edges[ti].vj;

#if 0
      printf("(%d %d) = %d\n", vi+1, vj+1, sups[ti]);
#endif

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        ei      = xaii[vi].end-1;
        eistart = xaii[vi].start;
        vik     = aii[ei].vj;

        ej      = xaii[vj].end-1;
        ejstart = xaii[vj].start;
        vjk     = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei >= eistart && ej >= ejstart) {
          if (vik > vjk) {
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;
          }
          else if (vjk > vik) {
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;
          }
          else {
            updindices[nupdates++] = ids[ei];
            updindices[nupdates++] = ids[ej];

            sups[ids[ei]]--;
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;

            sups[ids[ej]]--;
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    /* update the shead[k] information, for the subsequent updates */
    shead[k] = ti;
    slist[ti].peid = -1;


    /* add up sups[:] */
    int64_t total_sup = 0;
    #pragma omp parallel for schedule(static) reduction(+:total_sup)
    for(int64_t e = 0; e < nedges; ++e) {
      if(sups[e] >= 0) {
        total_sup += sups[e];
      }
    }
#if VERBOSE
    printf("  edges-left: %7zd (%5.2f%%), total-support: %7zd\n",
        nleft, 100. * (double)nleft / (double)nedges, total_sup);
#endif

    if (nupdates > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdates; ei++) {
        ti = updindices[ei];

        if (sups[ti] < 0 || sups[ti] == slist[ti].sup)
          continue; /* we have already deleted or updated this */

        /* remove ti from its current list */
        sup = (slist[ti].sup <= k ? k : slist[ti].sup);  /* see the comment in the "add" */
        if (shead[sup] != ti) { /* if ti was not the head */
          slist[slist[ti].peid].neid = slist[ti].neid;
          slist[slist[ti].neid].peid = slist[ti].peid;
        }
        else {
          shead[sup] = slist[ti].neid;
          slist[slist[ti].neid].peid = -1;
        }

        /* add ti to the head of the new list */
        sup = (sups[ti] <= k ? k : sups[ti]);  /* put all the <k support into the support
                                                  list that we are currently operating on */
        slist[ti].sup  = sups[ti];
        slist[ti].peid = -1;
        slist[ti].neid = shead[sup];
        slist[shead[sup]].peid = ti;
        shead[sup] = ti;
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }

    gk_stopwctimer(timer_currk);

    /* add up sups[:] */
    total_sup = 0;
    #pragma omp parallel for schedule(static) reduction(+:total_sup)
    for(int64_t e = 0; e < nedges; ++e) {
      if(sups[e] >= 0) {
        total_sup += sups[e];
      }
    }

    printf("k: %7d; edges-left: %7zd (%5.2f%%), total-support: %7zd, "
           "nltriangles: %7d, time (s): %6.3f\n",
        k+2, nleft, 100. * (double)nleft / (double)nedges, total_sup,
        nltriangles, timer_currk);
    gk_clearwctimer(timer_currk);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  slist--;
  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &shead, &slist, &updindices, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    list-traversal + bucket-based supports + O(1) deletion

    Uses an in-place support-based bucket implementation in the form of
    a linked list.

    This is an experiment in concurrency by pre-marking the edges that
    will be removed, followed by support updating, and finally actual
    edge removal.
*/
/*************************************************************************/
int64_t kt_Baseline5c(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  struct slist_s {
    ssize_t neid, peid;
    int32_t sup;
  } *slist;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nltriangles, sup;
  ssize_t ci, ti, ei, eistart, eiend, ej, ejstart, ejend, tsup0=0, tsup, ncand;
  int64_t nedges, nleft, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  int32_t k, nsups, *sups;
  ssize_t *ids, *shead;
  ssize_t nupdates, nmaxupdates, *updindices;
  gk_i64kv_t *cand;

  double timer_currk = 0.;

  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  /* the heads and "link list" that form the support buckets */
  shead = gk_zsmalloc(nsups, -1, "shead");
  slist = (struct slist_s *)gk_malloc((nedges+1)*sizeof(struct slist_s), "slist");

  slist++;  /* this is to allow slist[-1] to be valid */
  for (ei=0; ei<nedges; ei++) {
    slist[ei].sup  = sups[ei];
    slist[ei].peid = -1;
    slist[ei].neid = shead[sups[ei]];
    if (shead[sups[ei]] != -1)
      slist[shead[sups[ei]]].peid = ei;
    shead[sups[ei]] = ei;
  }

  nmaxupdates = nedges + 2*nvtxs;
  updindices = gk_zmalloc(nmaxupdates, "updindices");

  cand = gk_i64kvmalloc(nedges, "cand");

  gk_stopwctimer(vault->timer_ktsetup);

  if (params->dbglvl&1) {
    for (tsup0=0, ei=0; ei<nedges; ei++)
      tsup0 += (sups[ei] > 0 ? sups[ei] : 0);
  }

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;
  nleft      = nedges;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && nleft>0; k++) {
    nltriangles = 0;

    gk_clearwctimer(timer_currk);
    gk_startwctimer(timer_currk);
BACK:

    nupdates = 0;
    for (ncand=0, ti=shead[k]; ti!=-1; ti=slist[ti].neid, ncand++) {
      if (nupdates + 2*nvtxs > nmaxupdates)
        break;

      vi = edges[ti].vi;
      vj = edges[ti].vj;

      cand[ncand].key = vj*nvtxs + (nvtxs-vi);
      cand[ncand].val = ti;

      nupdates += 2*sups[ti];
    }

    /* update the shead[k] information, for the subsequent updates */
    shead[k] = ti;
    slist[ti].peid = -1;

    /* sort the edges to be removed in a "nice" order */
    gk_i64kvsorti(ncand, cand);

    nupdates = 0;
    for (ci=0; ci<ncand; ci++) {
      ti = cand[ci].val;

      nleft--;

      vi = edges[ti].vi;
      vj = edges[ti].vj;

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        ei      = xaii[vi].end-1;
        eistart = xaii[vi].start;
        vik     = aii[ei].vj;

        ej      = xaii[vj].end-1;
        ejstart = xaii[vj].start;
        vjk     = aii[ej].vj;

        /* decrease the support of the intersection */
        while (ei >= eistart && ej >= ejstart) {
          if (vik > vjk) {
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;
          }
          else if (vjk > vik) {
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;
          }
          else {
            updindices[nupdates++] = ids[ei];
            updindices[nupdates++] = ids[ej];

            sups[ids[ei]]--;
            ei  -= aii[ei].dec;
            vik  = aii[ei].vj;

            sups[ids[ej]]--;
            ej  -= aii[ej].dec;
            vjk  = aii[ej].vj;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (nupdates > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdates; ei++) {
        ti = updindices[ei];

        if (sups[ti] < 0 || sups[ti] == slist[ti].sup)
          continue; /* we have already deleted or updated this */

        /* remove ti from its current list */
        sup = (slist[ti].sup <= k ? k : slist[ti].sup);  /* see the comment in the "add" */
        if (shead[sup] != ti) { /* if ti was not the head */
          slist[slist[ti].peid].neid = slist[ti].neid;
          slist[slist[ti].neid].peid = slist[ti].peid;
        }
        else {
          shead[sup] = slist[ti].neid;
          slist[slist[ti].neid].peid = -1;
        }

        /* add ti to the head of the new list */
        sup = (sups[ti] <= k ? k : sups[ti]);  /* put all the <k support into the support
                                                  list that we are currently operating on */
        slist[ti].sup  = sups[ti];
        slist[ti].peid = -1;
        slist[ti].neid = shead[sup];
        slist[shead[sup]].peid = ti;
        shead[sup] = ti;
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }
    gk_stopwctimer(timer_currk);

    /* add up sups[:] */
    int64_t total_sup = 0;
    #pragma omp parallel for schedule(static) reduction(+:total_sup)
    for(int64_t e = 0; e < nedges; ++e) {
      if(sups[e] >= 0) {
        total_sup += sups[e];
      }
    }
    printf("k: %7d; edges-left: %7zd (%5.2f%%), total-support: %7zd, "
           "nltriangles: %7d, time (s): %6.3f\n",
        k+2, nleft, 100. * (double)nleft / (double)nedges, total_sup,
        nltriangles, timer_currk);


    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  slist--;
  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &shead, &slist, &updindices, &cand, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline6(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eip, eiend, ej, ejp, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  int32_t k, sk, nsups, *sups;
  ssize_t si, eid, *ids, *sptr, *sind, *sloc;
  ssize_t pi, hi, hiend, key, htsize, *htablek, *htablev;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        ids[xaii[vj].start]     = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        ids[xaii[vi].start]     = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  /* setup the support buckets and all associated information */
  gk_startwctimer(vault->timer_1);
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  sptr = gk_zsmalloc(nsups+1, 0, "sptr");
  sind = gk_zmalloc(nedges, "sind");
  sloc = gk_zmalloc(nedges, "sloc");

  for (ei=0; ei<nedges; ei++)
    sptr[sups[ei]]++;
  MAKECSR(vi, nsups, sptr);
  for (ei=0; ei<nedges; ei++) {
    sloc[ei] = sptr[sups[ei]];
    sind[sptr[sups[ei]]++] = ei;
  }
  SHIFTCSR(vi, nsups, sptr);
  gk_stopwctimer(vault->timer_1);

  /* create and populate the hash-table for the edges */
  htsize = 11*nedges;
  htablek = gk_zsmalloc(htsize, -1, "htablek");
  htablev = gk_zmalloc(htsize, "htablev");
  for (vi=0; vi<nvtxs; vi++) {
    for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei++) {
      vj  = aii[ei].vj;
      key = vi*nvtxs+vj;
      for (pi=0; pi<1000; pi++) {
        hi = hfun1(vi, vj, pi, htsize);
        if (htablek[hi] == -1) {
          htablek[hi] = key;
          htablev[hi] = ei;
          //printf("%d %d %zd, %zd %zd\n", vi, vj, pi, ei, hi);
          break;
        }
      }
      GKASSERT(pi < 1000);
    }
  }
  printf("htsize: %zd\n", htsize);

  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && sptr[k]!=nedges; k++) {
    nltriangles = 0;
    nlocal = 0;
    for (si=sptr[k]; si<sptr[k+1]; si++, nlocal++) {
      ti = sind[si];
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (sups[ti] > 0) {
        gk_startwctimer(vault->timer_3);
        sup = sups[ti];

        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, sups[ti],
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        /* use the htable to find intersections */
        for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei+=aii[ei].inc) {
          if ((vk = aii[ei].vj) == vj)
            continue;
          key = vj*nvtxs+vk;
          for (pi=0; pi<1000; pi++) {
            hi = hfun1(vj, vk, pi, htsize);
            if (htablek[hi] == -1 || htablek[hi] == key)
              break;
          }
          GKASSERT(pi<1000);

          if (htablek[hi] == key) {
            ej = htablev[hi];
            GKASSERT(aii[ei].vj == aii[ej].vj);

            if (sups[ids[ej]] < 0)
              continue;

            //printf("  ** vk: %d, ei: %zd, ej: %zd\n", vk, ei, ej);
            nltriangles++;

            eid = ids[ei];
            GKASSERT(sups[eid] >= 0);
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;

            eid = ids[ej];
            GKASSERT(sups[eid] >= 0);
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;

            if (--sup == 0)
              break;
          }
        }
        gk_stopwctimer(vault->timer_3);
      }

      /* remove the edge from both adjacency lists */
      gk_startwctimer(vault->timer_5);
      for (eip=-1, ei=xaii[vi].start; aii[ei].vj<vj; eip=ei, ei+=aii[ei].inc);
      if (eip == -1)
        xaii[vi].start += aii[ei].inc;
      else
        aii[eip].inc += aii[ei].inc;
      gk_stopwctimer(vault->timer_5);

      gk_startwctimer(vault->timer_6);
      for (ejp=-1, ej=xaii[vj].start; aii[ej].vj<vi; ejp=ej, ej+=aii[ej].inc);
      if (ejp == -1)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ejp].inc += aii[ej].inc;
      gk_stopwctimer(vault->timer_6);

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &sptr, &sind, &sloc, &htablek, &htablev, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline7(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vj, vi;
  } *edges;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  gk_ipq_t *queue;

  int32_t k, vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eiend, eistart, ej, ejend, nc, np;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt, *kadjncy;
  ssize_t nupdates, maxnupdates, *updindices;
  int32_t *sups;
  ssize_t eid, *ids;
  ssize_t *htsizes, *htablev;
  int32_t l, htsize, *htablek;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii    = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  kadjncy = gk_i32malloc(2*nedges+1, "kadjncy");
  edges   = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups    = gk_i32malloc(nedges, "sups");
  ids     = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        kadjncy[xaii[vj].start] = vi;
        ids[xaii[vj].start]     = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        kadjncy[xaii[vi].start] = adjncy[ei];
        ids[xaii[vi].start]     = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  /* create and populate the priority queue
     [the negative on the key is because my priority queues are max priority queues */
  gk_startwctimer(vault->timer_1);
  queue = gk_ipqCreate(nedges);
  for (ei=0; ei<nedges; ei++)
    gk_ipqInsert(queue, ei, -sups[ei]);

  maxnupdates = 2*nvtxs + nedges;
  updindices = gk_zmalloc(maxnupdates, "updindices");
  gk_stopwctimer(vault->timer_1);

  printf("Pqueue stats: nnodes: %zu, nedges: %zd\n", gk_ipqLength(queue), nedges);


  /* create and populate the hash-table for the edges */
  htsizes = gk_zmalloc(nvtxs+1, "htsizes");
  htsize = 0;
  for (vi=0; vi<nvtxs; vi++) {
    if ((xaii[vi].end-xaii[vi].start)<<2 > 1 + (htsize>>4) + (htsize>>1)) {
      htsize = xaii[vi].end-xaii[vi].start;
      for (l=1; htsize>(1<<l); l++);
      htsize = (1<<(l+3));
    }
    htsizes[vi] = htsize;
  }
  MAKECSR(vi, nvtxs, htsizes);
  printf("Size of htable: %zd (%5.2f)\n", htsizes[nvtxs], (float)htsizes[nvtxs]/xaii[nvtxs-1].end);

  htablek = gk_i32smalloc(htsizes[nvtxs], -1, "htablek");
  htablev = gk_zmalloc(htsizes[nvtxs], "htablev");
  for (nc=0, vi=0; vi<nvtxs; vi++) {
    htablek += htsizes[vi];
    htablev += htsizes[vi];
    htsize   = htsizes[vi+1]-htsizes[vi]-1;

    for (ei=xaii[vi].end-1, eistart=xaii[vi].start; ei>=eistart; ei--) {
      vk = kadjncy[ei];
      for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize, nc++);
      htablek[l] = vk;
      htablev[l] = ids[ei];
    }

    htablek -= htsizes[vi];
    htablev -= htsizes[vi];
  }

  printf("htable collision rate: %.4lf%%\n", (double)100.0*nc/xaii[nvtxs-1].end);


  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  nc = np = 0;
  for (k=1;;k++) {
    if (gk_ipqLength(queue) == 0)
      break;

    nltriangles = 0;
    nlocal = 0;

BACK:
    nupdates = 0;
    while (k >= -gk_ipqSeeTopKey(queue)) {
      if (nupdates + 2*nvtxs > maxnupdates)
        break;
      if (gk_ipqLength(queue) == 0)
        break;

      nlocal++;

      ti = gk_ipqGetTop(queue);
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (xaii[vi].end-xaii[vi].start > xaii[vj].end-xaii[vj].start) {
        vi = vj;
        vj = edges[ti].vi;
      }

      /* remove the edge from both adjacency lists */
      gk_startwctimer(vault->timer_5);
      for (ei=xaii[vi].start; kadjncy[ei]!=vj; ei++);
      GKASSERT(ei<xaii[vi].end);
      kadjncy[ei] = kadjncy[--xaii[vi].end];
      ids[ei]     = ids[xaii[vi].end];

      for (ej=xaii[vj].start; kadjncy[ej]!=vi; ej++);
      GKASSERT(ej<xaii[vj].end);
      kadjncy[ej] = kadjncy[--xaii[vj].end];
      ids[ej]     = ids[xaii[vj].end];
      gk_stopwctimer(vault->timer_5);

      if (sups[ti] > 0) {
        gk_startwctimer(vault->timer_3);
        sup = sups[ti];

        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, sups[ti],
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        /* use the htable to find intersections */
        htablek += htsizes[vj];
        htablev += htsizes[vj];
        htsize   = htsizes[vj+1]-htsizes[vj]-1;
        np += xaii[vi].end-xaii[vi].start;
        for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei++) {
          vk = kadjncy[ei];
          for (l=(vk&htsize); htablek[l]!=-1 && htablek[l]!=vk; l=(l+1)&htsize, nc++);
          if (htablek[l] == vk) {
            eid = htablev[l];

            /* the following is to deal with the case that htables contain all edges,
               even those that were already removed */
            if (sups[eid] < 0)
              continue;

            //printf("  ** vk: %d, ei: %zd, ej: %zd\n", vk, ei, ej);
            nltriangles++;

            sups[eid]--;
            updindices[nupdates++] = eid;

            eid = ids[ei];
            sups[eid]--;
            updindices[nupdates++] = eid;

            if (--sup == 0)
              break;
          }
        }
        htablek -= htsizes[vj];
        htablev -= htsizes[vj];
        gk_stopwctimer(vault->timer_3);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (nupdates > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdates; ei++) {
        if (sups[updindices[ei]] >= 0)
          if (queue->heap[queue->locator[updindices[ei]]].key != -sups[updindices[ei]])
            gk_ipqUpdate(queue, updindices[ei], -sups[updindices[ei]]);
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);
  printf("probe collision rate: %.4lf%% [np+nc: %zd]\n", (double)100.0*nc/np, np+nc);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_ipqDestroy(queue);

  gk_free((void **)&edges, &kadjncy, &xaii, &ids, &sups, &updindices,
      &htsizes, &htablek, &htablev, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline8(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vi, vj;
  } *edges;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eiend, eistart, ej, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt, *kadjncy;
  int32_t k, sk, nsups, *sups;
  ssize_t si, eid, *ids, *sptr, *sind, *sloc;
  ssize_t hiend, *htsizes, *htablev;
  int32_t l, htsize, *htablek;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii    = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  kadjncy = gk_i32malloc(2*nedges+1, "kadjncy");
  edges   = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups    = gk_i32malloc(nedges, "sups");
  ids     = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        kadjncy[xaii[vj].start] = vi;
        ids[xaii[vj].start]     = nedges;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        kadjncy[xaii[vi].start] = adjncy[ei];
        ids[xaii[vi].start]     = nedges;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  printf("Edges with non-zero support: %zd (%.2f%%)\n", nedges, (float)100.00*nedges/xadj[nvtxs]);

  /* setup the support buckets and all associated information */
  gk_startwctimer(vault->timer_1);
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  sptr = gk_zsmalloc(nsups+1, 0, "sptr");
  sind = gk_zmalloc(nedges, "sind");
  sloc = gk_zmalloc(nedges, "sloc");

  for (ei=0; ei<nedges; ei++)
    sptr[sups[ei]]++;
  MAKECSR(vi, nsups, sptr);
  for (ei=0; ei<nedges; ei++) {
    sloc[ei] = sptr[sups[ei]];
    sind[sptr[sups[ei]]++] = ei;
  }
  SHIFTCSR(vi, nsups, sptr);
  gk_stopwctimer(vault->timer_1);

  /* create and populate the hash-table for the edges */
  htsizes = gk_zmalloc(nvtxs+1, "htsizes");
  htsize = 0;
  for (vi=0; vi<nvtxs; vi++) {
    if ((xaii[vi].end-xaii[vi].start)<<2 > 1 + (htsize>>4) + (htsize>>1)) {
      htsize = xaii[vi].end-xaii[vi].start;
      for (l=1; htsize>(1<<l); l++);
      htsize = (1<<(l+3));
    }
    htsizes[vi] = htsize;
  }
  MAKECSR(vi, nvtxs, htsizes);
  printf("Size of htable: %zd (%5.2f)\n", htsizes[nvtxs], (float)htsizes[nvtxs]/xaii[nvtxs-1].end);

  htablek = gk_i32smalloc(htsizes[nvtxs], -1, "htablek");
  htablev = gk_zmalloc(htsizes[nvtxs], "htablev");
  for (vi=0; vi<nvtxs; vi++) {
    htablek += htsizes[vi];
    htablev += htsizes[vi];
    htsize   = htsizes[vi+1]-htsizes[vi]-1;

    for (ei=xaii[vi].end-1, eistart=xaii[vi].start; ei>=eistart; ei--) {
      vk = kadjncy[ei];
      for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize);
      htablek[l] = vk;
      htablev[l] = ids[ei];
    }

    htablek -= htsizes[vi];
    htablev -= htsizes[vi];
  }
  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && sptr[k]!=nedges; k++) {
    nltriangles = 0;
    nlocal = 0;
    for (si=sptr[k]; si<sptr[k+1]; si++, nlocal++) {
      ti = sind[si];
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (xaii[vi].end-xaii[vi].start > xaii[vj].end-xaii[vj].start) {
        vi = vj;
        vj = edges[ti].vi;
      }

      /* remove the edge from both adjacency lists */
      gk_startwctimer(vault->timer_5);
      for (ei=xaii[vi].start; kadjncy[ei]!=vj; ei++);
      GKASSERT(ei<xaii[vi].end);
      kadjncy[ei] = kadjncy[--xaii[vi].end];
      ids[ei]     = ids[xaii[vi].end];

      for (ej=xaii[vj].start; kadjncy[ej]!=vi; ej++);
      GKASSERT(ej<xaii[vj].end);
      kadjncy[ej] = kadjncy[--xaii[vj].end];
      ids[ej]     = ids[xaii[vj].end];
      gk_stopwctimer(vault->timer_5);

      if (sups[ti] > 0) {
        gk_startwctimer(vault->timer_3);
        sup = sups[ti];

        /*
        printf("%d %d %d %zd %zd %d\n",
            vi, vj, sups[ti],
            xaii[vi].end-xaii[vi].start,
            xaii[vj].end-xaii[vj].start,
            nltriangles);
        */

        /* use the htable to find intersections */
        htablek += htsizes[vj];
        htablev += htsizes[vj];
        htsize   = htsizes[vj+1]-htsizes[vj]-1;
        for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei++) {
          vk = kadjncy[ei];
          for (l=(vk&htsize); htablek[l]!=-1 && htablek[l]!=vk; l=(l+1)&htsize);
          if (htablek[l] == vk) {
            eid = htablev[l];

            /* the following is to deal with the case that htables contain all edges,
               even those that were already removed */
            if (sups[eid] < 0)
              continue;

            //printf("  ** vk: %d, ei: %zd, ej: %zd\n", vk, ei, ej);
            nltriangles++;

            //GKASSERT(sups[eid] >= 0);
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;


            eid = ids[ei];
            //GKASSERT(sups[eid] >= 0);
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;

            if (--sup == 0)
              break;
          }
        }
        htablek -= htsizes[vj];
        htablev -= htsizes[vj];
        gk_stopwctimer(vault->timer_3);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_free((void **)&edges, &kadjncy, &xaii, &ids, &sups, &sptr, &sind, &sloc,
      &htsizes, &htablek, &htablev, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.

*/
/*************************************************************************/
int64_t kt_Baseline10(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vj, vi;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eip, eiend, eistart, ej, ejp, ejend;
  int64_t nedges, ntriangles;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  int32_t k, sk, nsups, *sups;
  ssize_t si, eid, *ids, *sptr, *sind, *sloc;
  ssize_t hiend, *htsizes, *htablev;
  int32_t l, htsize, *htablek;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  sptr = gk_zsmalloc(nsups+1, 0, "sptr");
  sind = gk_zmalloc(nedges, "sind");
  sloc = gk_zmalloc(nedges, "sloc");

  for (ei=0; ei<nedges; ei++)
    sptr[sups[ei]]++;
  MAKECSR(vi, nsups, sptr);
  for (ei=0; ei<nedges; ei++) {
    sloc[ei] = sptr[sups[ei]];
    sind[sptr[sups[ei]]++] = ei;
  }
  SHIFTCSR(vi, nsups, sptr);

  /* create and populate the hash-table for the edges */
  htsizes = gk_zmalloc(nvtxs+1, "htsizes");
  htsize = 0;
  for (vi=0; vi<nvtxs; vi++) {
    if ((xaii[vi].end-xaii[vi].start)<<2 > 1 + (htsize>>4) + (htsize>>1)) {
      htsize = xaii[vi].end-xaii[vi].start;
      for (l=1; htsize>(1<<l); l++);
      htsize = (1<<(l+3));
    }
    htsizes[vi] = htsize;
  }
  MAKECSR(vi, nvtxs, htsizes);
  printf("Size of htable: %zd (%5.2f)\n", htsizes[nvtxs], (float)htsizes[nvtxs]/xaii[nvtxs-1].end);

  htablek = gk_i32smalloc(htsizes[nvtxs], -1, "htablek");
  htablev = gk_zmalloc(htsizes[nvtxs], "htablev");
  for (vi=0; vi<nvtxs; vi++) {
    htablek += htsizes[vi];
    htablev += htsizes[vi];
    htsize   = htsizes[vi+1]-htsizes[vi]-1;

    for (ei=xaii[vi].end-1, eistart=xaii[vi].start; ei>=eistart; ei--) {
      vk = aii[ei].vj;
      for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize);
      htablek[l] = vk;
      htablev[l] = ids[ei];
    }

    htablek -= htsizes[vi];
    htablev -= htsizes[vi];
  }
  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && sptr[k]!=nedges; k++) {
    nltriangles = 0;
    nlocal = 0;
    for (si=sptr[k]; si<sptr[k+1]; si++, nlocal++) {
      ti = sind[si];
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        /* use the htable to find intersections */
        htablek += htsizes[vj];
        htablev += htsizes[vj];
        htsize   = htsizes[vj+1]-htsizes[vj]-1;
        for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei+=aii[ei].inc) {
          vk = aii[ei].vj;
          for (l=(vk&htsize); htablek[l]!=-1 && htablek[l]!=vk; l=(l+1)&htsize);
          if (htablek[l] == vk) {
            eid = htablev[l];

            /* the following is to deal with the case that htables contain all edges,
               even those that were already removed */
            if (sups[eid] < 0)
              continue;

            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;


            eid = ids[ei];
            /* NOTE: if sups[eid] <= k, we do not move it, as it is in
               the bucket that we are currently processing */
            if ((sk=sups[eid]) > k) {
              /* put in eid's spot the first edge of that bucket */
              sind[sloc[eid]] = sind[sptr[sk]];
              /* update the sloc of the moved edge */
              sloc[sind[sloc[eid]]] = sloc[eid];
              /* put eid into the place from which we moved an edge */
              sind[sptr[sk]] = eid;
              /* update the sloc of eid */
              sloc[eid] = sptr[sk];
              /* advance the sptr of that bucket to reflect that it shrunk */
              sptr[sk]++;
            }
            sups[eid]--;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);
        htablek -= htsizes[vj];
        htablev -= htsizes[vj];
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    if (params->dbglvl&1 && nlocal > 0)
      printf("k: %7d; nlocal: %7d, nltriangles: %7d\n", k, nlocal, nltriangles);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &sptr, &sind, &sloc,
      &htsizes, &htablek, &htablev, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    hash-tables + bucket-based supports + O(1) deletion

    Uses an in-place support-based bucket implementation in the form of
    a linked list.
*/
/*************************************************************************/
int64_t kt_Baseline10b(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vj, vi;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
  } *xaii;

  struct slist_s {
    ssize_t neid, peid;
    int32_t sup;
  } *slist;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup;
  ssize_t ti, ei, eip, eiend, eistart, ej, ejp, ejend;
  int64_t nedges, ntriangles, nleft;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  ssize_t hiend, *htsizes, *htablev;
  int32_t l, htsize, *htablek;
  int32_t k, nsups, *sups;
  ssize_t eid, *ids, *shead;
  ssize_t nupdates, nmaxupdates, *updindices;

  double timer_currk = 0.;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].end = xaii[vi+1].start;

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  /* the heads and "link list" that form the support buckets */
  shead = gk_zsmalloc(nsups, -1, "shead");
  slist = (struct slist_s *)gk_malloc((nedges+1)*sizeof(struct slist_s), "slist");

  slist++;  /* this is to allow slist[-1] to be valid */
  for (ei=0; ei<nedges; ei++) {
    slist[ei].sup  = sups[ei];
    slist[ei].peid = -1;
    slist[ei].neid = shead[sups[ei]];
    if (shead[sups[ei]] != -1)
      slist[shead[sups[ei]]].peid = ei;
    shead[sups[ei]] = ei;
  }

  nmaxupdates = nedges + 2*nvtxs;
  updindices = gk_zmalloc(nmaxupdates, "updindices");


  /* create and populate the hash-table for the edges */
  htsizes = gk_zmalloc(nvtxs+1, "htsizes");
  htsize = 0;
  for (vi=0; vi<nvtxs; vi++) {
    if ((xaii[vi].end-xaii[vi].start)<<2 > 1 + (htsize>>4) + (htsize>>1)) {
      htsize = xaii[vi].end-xaii[vi].start;
      for (l=1; htsize>(1<<l); l++);
      htsize = (1<<(l+3));
    }
    htsizes[vi] = htsize;
  }
  MAKECSR(vi, nvtxs, htsizes);
  printf("Size of htable: %zd (%5.2f)\n", htsizes[nvtxs], (float)htsizes[nvtxs]/xaii[nvtxs-1].end);

  htablek = gk_i32smalloc(htsizes[nvtxs], -1, "htablek");
  htablev = gk_zmalloc(htsizes[nvtxs], "htablev");
  for (vi=0; vi<nvtxs; vi++) {
    htablek += htsizes[vi];
    htablev += htsizes[vi];
    htsize   = htsizes[vi+1]-htsizes[vi]-1;

    for (ei=xaii[vi].end-1, eistart=xaii[vi].start; ei>=eistart; ei--) {
      vk = aii[ei].vj;
      for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize);
      htablek[l] = vk;
      htablev[l] = ids[ei];
    }

    htablek -= htsizes[vi];
    htablev -= htsizes[vi];
  }
  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;
  nleft = nedges;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && nleft>0; k++) {
    nltriangles = 0;
    gk_clearwctimer(timer_currk);
    gk_startwctimer(timer_currk);

BACK:
    nupdates = 0;
    for (ti=shead[k]; ti!=-1; ti=slist[ti].neid) {
      if (nupdates + 2*nvtxs > nmaxupdates)
        break;

      nleft--;

      vi = edges[ti].vi;
      vj = edges[ti].vj;

      /* remove the edge from both adjacency lists */
      ei = edges[ti].eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      ej = edges[ti].eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        /* use the htable to find intersections */
        htablek += htsizes[vj];
        htablev += htsizes[vj];
        htsize   = htsizes[vj+1]-htsizes[vj]-1;
        for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei+=aii[ei].inc) {
          vk = aii[ei].vj;
          for (l=(vk&htsize); htablek[l]!=-1 && htablek[l]!=vk; l=(l+1)&htsize);
          if (htablek[l] == vk) {
            eid = htablev[l];

            /* the following is to deal with the case that htables contain all edges,
               even those that were already removed */
            if (sups[eid] < 0)
              continue;

            sups[eid]--;
            updindices[nupdates++] = eid;

            eid = ids[ei];
            sups[eid]--;
            updindices[nupdates++] = eid;

            if (--sup == 0)
              break;
          }
        }
        GKASSERT(sup == 0);

        htablek -= htsizes[vj];
        htablev -= htsizes[vj];
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }

    /* update the shead[k] information, for the subsequent updates */
    shead[k] = ti;
    slist[ti].peid = -1;

    if (nupdates > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdates; ei++) {
        ti = updindices[ei];

        if (sups[ti] < 0 || sups[ti] == slist[ti].sup)
          continue; /* we have already deleted or updated this */

        /* remove ti from its current list */
        sup = (slist[ti].sup <= k ? k : slist[ti].sup);  /* see the comment in the "add" */
        if (shead[sup] != ti) { /* if ti was not the head */
          slist[slist[ti].peid].neid = slist[ti].neid;
          slist[slist[ti].neid].peid = slist[ti].peid;
        }
        else {
          shead[sup] = slist[ti].neid;
          slist[slist[ti].neid].peid = -1;
        }

        /* add ti to the head of the new list */
        sup = (sups[ti] <= k ? k : sups[ti]);  /* put all the <k support into the support
                                                  list that we are currently operating on */
        slist[ti].sup  = sups[ti];
        slist[ti].peid = -1;
        slist[ti].neid = shead[sup];
        slist[shead[sup]].peid = ti;
        shead[sup] = ti;
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }

    gk_stopwctimer(timer_currk);

    /* add up sups[:] */
    int64_t total_sup = 0;
    #pragma omp parallel for schedule(static) reduction(+:total_sup)
    for(int64_t e = 0; e < nedges; ++e) {
      if(sups[e] >= 0) {
        total_sup += sups[e];
      }
    }

    printf("k: %7d; edges-left: %7zd (%5.2f%%), total-support: %7zd, "
           "nltriangles: %7d, time (s): %6.3f\n",
        k+2, nleft, 100. * (double)nleft / (double)nedges, total_sup,
        nltriangles, timer_currk);


    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  slist--;
  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &slist, &shead, &updindices,
      &htsizes, &htablek, &htablev, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! This is the initial baseline version of k-truss decomposition.
    hash-tables + bucket-based supports + O(1) deletion

    Uses an in-place support-based bucket implementation in the form of
    a linked list.
*/
/*************************************************************************/
int64_t kt_Baseline10g(params_t *params, vault_t *vault)
{
  struct edge_s {
    int32_t vj, vi;
    ssize_t eij, eji;
  } *edges;

  struct aii_s {
    int32_t vj;
    int32_t inc, dec;
  } *aii;

  struct xaii_s {
    int64_t start, end;
    int32_t degree, lastdegree;
    int32_t htsize, nc;
    ssize_t *htablev;
    int32_t *htablek;
  } *xaii;

  struct slist_s {
    ssize_t neid, peid;
    int32_t sup;
  } *slist;

  int32_t vi, vik, vj, vjk, vk, nvtxs, nlocal, nltriangles, sup, nc;
  ssize_t ci, ti, ei, eiend, eistart, ej, ejstart, ejend, htsize, eij, eji;
  int64_t nedges, ntriangles, nleft, ncand;
  ssize_t *xadj;
  int32_t *adjncy, *adjwgt;
  ssize_t hiend, *htablev;
  int32_t l, *htablek;
  int32_t k, nsups, *sups;
  ssize_t eid, *ids, *shead;
  ssize_t nupdhead, nupdtail, *updindices, *updmarker;
  gk_i64kv_t *cand;


  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  adjwgt = vault->ugraph->iadjwgt = gk_i32smalloc(xadj[nvtxs], 0, "adjwgt");
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);


  gk_startwctimer(vault->timer_ktsetup);
  /* determine the number of edges with non-zero support */
  for (nedges=0, ei=0, eiend=xadj[nvtxs]; ei<eiend; ei++) {
    if (adjwgt[ei] > 0)
      nedges++;
  }

  /* allocate memory for the adjancency lists, which in addition to the
     adjancent vertex it will store the decrement (for skip-list) and
     the ID for priority queue */
  xaii  = (struct xaii_s *)gk_malloc((nvtxs+1)*sizeof(struct xaii_s), "xaii");
  aii   = (struct aii_s *)gk_malloc((2*nedges+1)*sizeof(struct aii_s), "aii");
  edges = (struct edge_s *)gk_malloc((nedges+1)*sizeof(struct edge_s), "edges");
  sups  = gk_i32malloc(nedges, "sups");
  ids   = gk_zmalloc(2*nedges+1, "ids");

  for (vi=0; vi<nvtxs; vi++)
    xaii[vi].start = 0;

  /* determine sizes */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        xaii[vi].start++;
        xaii[adjncy[ei]].start++;

        edges[nedges].vi = vi;
        edges[nedges].vj = adjncy[ei];
        sups[nedges]     = adjwgt[ei];
        nedges++;
      }
    }
  }
  /* the MAKECSR equivalent */
  for (vi=1; vi<nvtxs; vi++)
    xaii[vi].start += xaii[vi-1].start;
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* populate it into two steps to ensure that the sorted order is maintained */
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        vj = adjncy[ei];
        aii[xaii[vj].start].vj  = vi;
        aii[xaii[vj].start].inc = 1;
        aii[xaii[vj].start].dec = 1;
        ids[xaii[vj].start]     = nedges;
        edges[nedges].eji       = xaii[vj].start;
        xaii[vj].start++;
        nedges++;
      }
    }
  }
  for (nedges=0, vi=0; vi<nvtxs; vi++) {
    for (ei=xadj[vi], eiend=xadj[vi+1]; ei<eiend; ei++) {
      if (adjwgt[ei] > 0) {
        aii[xaii[vi].start].vj  = adjncy[ei];
        aii[xaii[vi].start].inc = 1;
        aii[xaii[vi].start].dec = 1;
        ids[xaii[vi].start]     = nedges;
        edges[nedges].eij       = xaii[vi].start;
        xaii[vi].start++;
        nedges++;
      }
    }
  }
  /* the SHIFTCSR equivalent */
  for (vi=nvtxs; vi>0; vi--)
    xaii[vi].start = xaii[vi-1].start;
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for (vi=0; vi<nvtxs; vi++) {
    xaii[vi].end = xaii[vi+1].start;
    xaii[vi].degree = xaii[vi].lastdegree = xaii[vi].end-xaii[vi].start;
  }

#ifdef XXX
  /* free vault->ugraph */
  /* NOTE: I need to re-do how the decomposition is extracted in order to be able to
   * do this */
  gk_graph_Free(&vault->ugraph);
#endif

  /* setup the support buckets and all associated information */
  nsups = gk_i32max(nedges, sups, 1) + 1;
  printf("nsups: %d\n", nsups);

  /* the heads and "link list" that form the support buckets */
  shead = gk_zsmalloc(nsups, -1, "shead");
  slist = (struct slist_s *)gk_malloc((nedges+1)*sizeof(struct slist_s), "slist");

  slist++;  /* this is to allow slist[-1] to be valid */
  for (ei=0; ei<nedges; ei++) {
    slist[ei].sup  = sups[ei];
    slist[ei].peid = -1;
    slist[ei].neid = shead[sups[ei]];
    slist[shead[sups[ei]]].peid = ei;
    shead[sups[ei]] = ei;
  }

  updindices = gk_zmalloc(nedges, "updindices");
  updmarker  = gk_zsmalloc(nedges, -1, "updmarker");

  /* create and populate the hash-table for the edges */
  for (l=1, vi=0; vi<nvtxs; vi++) {
    htsize = xaii[vi].degree + (xaii[vi].degree>>1) + 1;
    for (; htsize>(1<<l); l++);
    xaii[vi].htsize = (1<<(l+0));
  }
  for (htsize=0, vi=0; vi<nvtxs; vi++)
    htsize += xaii[vi].htsize;
  printf("Size of htable: %zd (%5.2f)\n", htsize, (float)htsize/xaii[nvtxs-1].end);

  xaii[0].htablek = gk_i32smalloc(htsize, -1, "htablek");
  xaii[0].htablev = gk_zmalloc(htsize, "htablev");
  for (vi=1; vi<nvtxs; vi++) {
    xaii[vi].htablek = xaii[vi-1].htablek + xaii[vi-1].htsize;
    xaii[vi].htablev = xaii[vi-1].htablev + xaii[vi-1].htsize;
  }

  for (vi=0; vi<nvtxs; vi++) {
    htablek = xaii[vi].htablek;
    htablev = xaii[vi].htablev;
    htsize  = xaii[vi].htsize-1;

    for (nc=0, ei=xaii[vi].end-1, eistart=xaii[vi].start; ei>=eistart; ei--) {
      vk = aii[ei].vj;
      for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize, nc++);
      htablek[l] = vk;
      htablev[l] = ids[ei];
    }

    xaii[vi].nc = nc;
  }

  cand = gk_i64kvmalloc(nedges, "cand");

  gk_stopwctimer(vault->timer_ktsetup);

  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;
  nleft = nedges;

  gk_startwctimer(vault->timer_ktpeeling);
  /* get into the k-truss enumeration loop */
  for (k=1; k<nsups && nleft>0; k++) {
    nltriangles = 0;
    nupdtail = nedges;

BACK:
    gk_startwctimer(vault->timer_3);
    for (ncand=0, ti=shead[k]; ti!=-1; ti=slist[ti].neid, ncand++) {
      vi = edges[ti].vi;
      vj = edges[ti].vj;

      if (xaii[vi].degree > xaii[vj].degree) {
        vi = edges[ti].vj;
        vj = edges[ti].vi;
      }

      cand[ncand].key = (int64_t)vj*nvtxs + (int64_t)(nvtxs-vi);
      //cand[ncand].key = (int64_t)vj*nvtxs + (int64_t)vi;
      cand[ncand].val = (vi == edges[ti].vi ? ti : -ti);
    }

    /* update the shead[k] information, for the subsequent updates */
    shead[k] = ti;  // this should be -1
    slist[ti].peid = -1;

    /* sort the edges to be removed in a "nice" order */
    gk_i64kvsorti(ncand, cand);

    nupdhead = 0;
    for (ci=0; ci<ncand; ci++) {
      ti = cand[ci].val;

      nleft--;

      if (ti > 0) {  /* the sign determines the order */
        vi  = edges[ti].vi;
        vj  = edges[ti].vj;
        eij = edges[ti].eij;
        eji = edges[ti].eji;
      }
      else {
        ti  = -ti;
        vi  = edges[ti].vj;
        vj  = edges[ti].vi;
        eij = edges[ti].eji;
        eji = edges[ti].eij;
      }

      /* remove the edge from both adjacency lists */
      xaii[vi].degree--;
      //ei = edges[ti].eij;
      ei = eij;
      if (ei == xaii[vi].start)
        xaii[vi].start += aii[ei].inc;
      else
        aii[ei-aii[ei].dec].inc += aii[ei].inc;
      if (ei == xaii[vi].end-1)
        xaii[vi].end -= aii[ei].dec;
      else
        aii[ei+aii[ei].inc].dec += aii[ei].dec;

      xaii[vj].degree--;
      //ej = edges[ti].eji;
      ej = eji;
      if (ej == xaii[vj].start)
        xaii[vj].start += aii[ej].inc;
      else
        aii[ej-aii[ej].dec].inc += aii[ej].inc;
      if (ej == xaii[vj].end-1)
        xaii[vj].end -= aii[ej].dec;
      else
        aii[ej+aii[ej].inc].dec += aii[ej].dec;


      /* see if there is a need to rehash vj */
      if (xaii[vj].htsize >= 2048 && xaii[vj].nc && xaii[vj].degree < .75*xaii[vj].lastdegree) {
        htsize = xaii[vj].degree + (xaii[vj].degree>>1) + 1;
        for (l=1; htsize>(1<<l); l++);
        xaii[vj].htsize = gk_min(xaii[vj].htsize, (1<<(l+3)));
        xaii[vj].lastdegree = xaii[vj].degree;

        htablek = xaii[vj].htablek;
        htablev = xaii[vj].htablev;
        htsize  = xaii[vj].htsize-1;

        gk_i32set(htsize+1, -1, htablek); /* we only need to reset htablek */
        for (nc=0, ej=xaii[vj].end-1, ejstart=xaii[vj].start; ej>=ejstart; ej-=aii[ej].dec) {
          vk = aii[ej].vj;
          for (l=(vk&htsize); htablek[l]!=-1; l=(l+1)&htsize, nc++);
          htablek[l] = vk;
          htablev[l] = ids[ej];
        }

        xaii[vj].nc = nc;
      }


      if (sups[ti] > 0) {
        sup = sups[ti];

        nltriangles += sup;

        /* use the htable to find intersections */
        htablek = xaii[vj].htablek;
        htablev = xaii[vj].htablev;
        htsize  = xaii[vj].htsize-1;
        if (xaii[vj].nc) {
          for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei+=aii[ei].inc) {
            vk = aii[ei].vj;
            for (l=(vk&htsize); htablek[l]!=-1 && htablek[l]!=vk; l=(l+1)&htsize);
            if (htablek[l] == vk) {
              eid = htablev[l];

              /* the following is to deal with the case that htables contain all edges,
                 even those that were already removed */
              if (sups[eid] < 0)
                continue;

              sups[eid]--;
              if (sups[eid] <= k) { /* goes to the head */
                if (updmarker[eid] == -1) {
                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
                else if (updmarker[eid] > nupdhead) { /* move it from the tail */
                  updindices[updmarker[eid]]      = updindices[nupdtail];
                  updmarker[updindices[nupdtail]] = updmarker[eid];
                  nupdtail++;

                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
              }
              else { /* goes to the tail */
                if (updmarker[eid] == -1) {
                  updmarker[eid] = --nupdtail;
                  updindices[nupdtail] = eid;
                }
              }

              eid = ids[ei];
              sups[eid]--;
              if (sups[eid] <= k) { /* goes to the head */
                if (updmarker[eid] == -1) {
                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
                else if (updmarker[eid] > nupdhead) { /* move it from the tail */
                  updindices[updmarker[eid]]      = updindices[nupdtail];
                  updmarker[updindices[nupdtail]] = updmarker[eid];
                  nupdtail++;

                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
              }
              else { /* goes to the tail */
                if (updmarker[eid] == -1) {
                  updmarker[eid] = --nupdtail;
                  updindices[nupdtail] = eid;
                }
              }

              if (--sup == 0)
                break;
            }
          }
        }
        else { /* no collisions */
          for (ei=xaii[vi].start, eiend=xaii[vi].end; ei<eiend; ei+=aii[ei].inc) {
            vk = aii[ei].vj;
            if (htablek[vk&htsize] == vk) {
              eid = htablev[vk&htsize];

              /* the following is to deal with the case that htables contain all edges,
                 even those that were already removed */
              if (sups[eid] < 0)
                continue;

              sups[eid]--;
              if (sups[eid] <= k) { /* goes to the head */
                if (updmarker[eid] == -1) {
                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
                else if (updmarker[eid] > nupdhead) { /* move it from the tail */
                  updindices[updmarker[eid]]      = updindices[nupdtail];
                  updmarker[updindices[nupdtail]] = updmarker[eid];
                  nupdtail++;

                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
              }
              else { /* goes to the tail */
                if (updmarker[eid] == -1) {
                  updmarker[eid] = --nupdtail;
                  updindices[nupdtail] = eid;
                }
              }

              eid = ids[ei];
              sups[eid]--;
              if (sups[eid] <= k) { /* goes to the head */
                if (updmarker[eid] == -1) {
                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
                else if (updmarker[eid] > nupdhead) { /* move it from the tail */
                  updindices[updmarker[eid]]      = updindices[nupdtail];
                  updmarker[updindices[nupdtail]] = updmarker[eid];
                  nupdtail++;

                  updmarker[eid]         = nupdhead;
                  updindices[nupdhead++] = eid;
                }
              }
              else { /* goes to the tail */
                if (updmarker[eid] == -1) {
                  updmarker[eid] = --nupdtail;
                  updindices[nupdtail] = eid;
                }
              }

              if (--sup == 0)
                break;
            }
          }
        }
        GKASSERT(sup == 0);
      }

      sups[ti] = -k;  /* this is used for encoding the maximal value of k of that edge */
    }
    gk_stopwctimer(vault->timer_3);

    if (nupdhead > 0) {
      gk_startwctimer(vault->timer_4);
      for (ei=0; ei<nupdhead; ei++) {
        ti = updindices[ei];

        if (sups[ti] < 0)
          continue; /* we have already deleted this */

        updmarker[ti] = -1;

        ASSERT(sups[ti] <= k);
        ASSERT(slist[ti].sup > k);

        /* remove ti from its current list */
        sup = slist[ti].sup;
        if (shead[sup] != ti) { /* if ti was not the head */
          slist[slist[ti].peid].neid = slist[ti].neid;
          slist[slist[ti].neid].peid = slist[ti].peid;
        }
        else {
          shead[sup] = slist[ti].neid;
          slist[slist[ti].neid].peid = -1;
        }

        /* add ti to the head of the new list */
        slist[ti].sup  = sups[ti];
        slist[ti].peid = -1;
        slist[ti].neid = shead[k];
        slist[shead[k]].peid = ti;
        shead[k] = ti;
      }
      gk_stopwctimer(vault->timer_4);
      goto BACK;
    }

    if (nupdtail < nedges) {
      gk_startwctimer(vault->timer_5);
      for (ei=nupdtail; ei<nedges; ei++) {
        ti = updindices[ei];
        updmarker[ti] = -1;

        ASSERT(sups[ti] > k);

        /* remove ti from its current list */
        sup = slist[ti].sup;
        if (shead[sup] != ti) { /* if ti was not the head */
          slist[slist[ti].peid].neid = slist[ti].neid;
          slist[slist[ti].neid].peid = slist[ti].peid;
        }
        else {
          shead[sup] = slist[ti].neid;
          slist[slist[ti].neid].peid = -1;
        }

        /* add ti to the head of the new list */
        sup = sups[ti];
        slist[ti].sup  = sup;
        slist[ti].peid = -1;
        slist[ti].neid = shead[sup];
        slist[shead[sup]].peid = ti;
        shead[sup] = ti;
      }
      gk_stopwctimer(vault->timer_5);
    }

    if (params->dbglvl&1)
      printf("k: %7d; nleft: %7zd, nltriangles: %9d, nupd/nleft: %5.2f%%\n",
          k, nleft, nltriangles, 100.0*(nedges-nupdtail)/nleft);

    ntriangles += nltriangles;
  }
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("#triangles after peeling: %zd\n", ntriangles);

  /* create the output of the decomposition */
  kt_Sups2KTEdges(params, vault, k-1, sups);

  slist--;
  htablek = xaii[0].htablek;
  htablev = xaii[0].htablev;
  gk_free((void **)&edges, &aii, &xaii, &ids, &sups, &slist, &shead, &updindices,
      &updmarker, &htablek, &htablev, &cand, LTERM);

  return ntriangles;
}


/*************************************************************************/
/*! The hash-map-based edge-triangle-support counting routine that uses
    the JIK triangle enumeration scheme.

    This is the mapjikv2 tc version.
*/
/*************************************************************************/
int64_t kt_ComputeEdgeSupportPar(params_t *params, vault_t *vault)
{
  int32_t vi, vj, vk, vl, nvtxs, nlocal;
  ssize_t ei, eiend, ej, ejstart, ejend;
  int64_t ntriangles, ntriangles2;
  ssize_t *xadj, *txadj;
  int32_t *adjncy, *tadjncy, *adjwgt;
  int32_t l, tnc, nc, hmsize, tlsize, tlstart;

  gk_startwctimer(vault->timer_2);

  nvtxs  = vault->ugraph->nvtxs;
  xadj   = vault->ugraph->xadj;
  adjncy = vault->ugraph->adjncy;
  adjwgt = vault->ugraph->iadjwgt;

  txadj   = vault->lgraph->xadj;
  tadjncy = vault->lgraph->adjncy;

  /* determine the size of the hash-map and convert it into a format
     that is compatible with a bitwise AND operation */
  for (hmsize=0, vi=0; vi<nvtxs; vi++)
    hmsize = gk_max(hmsize, (int32_t)(xadj[vi+1]-xadj[vi]));
  for (l=1; hmsize>(1<<l); l++);
  hmsize = (1<<(l+4))-1;
  printf("& compatible maximum hmsize: %"PRId32"\n", hmsize);

  /* determine the size of the tail-map and allocate memory for it */
  for (vi=(nvtxs>>2); vi<nvtxs; vi++) {
    if ((txadj[vi+1]-txadj[vi])<<9 > vi)
      break;
    if ((xadj[vi+1]-xadj[vi])<<4 > nvtxs-vi)
      break;
  }
  tlsize  = nvtxs - vi + 100;
  tlstart = nvtxs-tlsize;
  printf("tlsize: %"PRId32"\n", tlsize);
  printf("tlstart: %"PRId32"\n", tlstart);

  /* start counting triangles */
  if (params->dbglvl&1)
    gk_startwctimer(vault->timer_4);

  /* use a combination of hmap and tmap */
  ntriangles = 0;
  tnc        = 0;
#pragma omp parallel default(none) shared(xadj, txadj, hmsize, params, tlstart, tlsize, adjncy, tadjncy, adjwgt) private(nc, ej, ejstart, ejend, l, nlocal, vi, vk, ei, eiend) reduction(+:ntriangles) reduction(+:tnc)
  {

    int32_t *hmap = gk_i32smalloc(hmsize+1, -1, "hmap");
    int32_t *tmap = gk_i32smalloc(tlsize, -1, "tmap");
    tmap -= tlstart; /* make indexing simpler */
    int32_t hmsizel = 0;

    #pragma omp for schedule(dynamic, DYNAMIC_CHUNK)
    for (vj=1; vj<tlstart; vj++) {
      if (xadj[vj+1] == xadj[vj] || txadj[vj+1] == txadj[vj])
        continue;

      /* if needed, increase the working hmsize */
      if ((xadj[vj+1]-xadj[vj])<<3 > 1 + (hmsizel>>4) + (hmsizel>>1)) {
        hmsizel = xadj[vj+1]-xadj[vj];
        for (l=1; hmsizel>(1<<l); l++);
        hmsizel = (1<<(l+4))-1;
      }

      /* hash Adj(vj) using hmap for the front and tmap for the last tlsize indices */
      for (nc=0, ej=ejstart=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
        if ((vk = adjncy[ej]) >= tlstart)
          break;
        for (l=(vk&hmsizel); hmap[l]!=-1; l=((l+1)&hmsizel), nc++);
        hmap[l] = ej-ejstart;
      }
      for (; ej<ejend; ej++)
        tmap[adjncy[ej]] = ej-ejstart;

      /* find intersections */
      if (nc > 0) { /* we had collisions */
        for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
          vi = tadjncy[ej];
          for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
            if ((vk = adjncy[ei]) >= tlstart)
              break;
            l = vk&hmsizel;
            if (hmap[l] == -1)
              continue;
            if (adjncy[ejstart+hmap[l]] == vk) {
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+hmap[l]]++;
              nlocal++;
              continue;
            }
            for (l=((l+1)&hmsizel); hmap[l]!=-1 && adjncy[ejstart+hmap[l]]!=vk; l=((l+1)&hmsizel));
            if (hmap[l]!=-1 && adjncy[ejstart+hmap[l]] == vk) {
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+hmap[l]]++;
              nlocal++;
            }
          }
          for (; ei<eiend; ei++) {
            if (tmap[adjncy[ei]] != -1) {
              assert(adjncy[ejstart+tmap[adjncy[ei]]] == adjncy[ei]);
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+tmap[adjncy[ei]]]++;
              nlocal++;
            }
          }

          if (nlocal > 0) {
            ntriangles += nlocal;

            assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
            #pragma omp atomic
            adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
          }
        }

        /* reset hmap/tmap */
        for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
          if ((vk = adjncy[ej]) >= tlstart)
            break;
          for (l=(vk&hmsizel); hmap[l]==-1 || adjncy[ejstart+hmap[l]]!=vk; l=((l+1)&hmsizel));
          hmap[l] = -1;
        }
        for (; ej<ejend; ej++)
          tmap[adjncy[ej]] = -1;
      }
      else { /* there were no collisons */
        for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
          vi = tadjncy[ej];

          for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
            if ((vk = adjncy[ei]) >= tlstart)
              break;
            if (hmap[vk&hmsizel]!=-1 && adjncy[ejstart+hmap[vk&hmsizel]] == vk) {
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+hmap[vk&hmsizel]]++;
              nlocal++;
            }
          }
          for (; ei<eiend; ei++) {
            if (tmap[adjncy[ei]] != -1) {
              assert(adjncy[ejstart+tmap[adjncy[ei]]] == adjncy[ei]);
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+tmap[adjncy[ei]]]++;
              nlocal++;
            }
          }

          if (nlocal > 0) {
            ntriangles += nlocal;

            assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
            #pragma omp atomic
            adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
          }
        }

        /* reset hmap/tmap */
        for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++) {
          if ((vk = adjncy[ej]) >= tlstart)
            break;
          hmap[vk&hmsizel] = -1;
        }
        for (; ej<ejend; ej++)
          tmap[adjncy[ej]] = -1;
      }
    }

    tmap += tlstart;
    gk_free((void **)&hmap, &tmap, LTERM);
  }

  int32_t tlstart_idx = tlstart;
  if (tlstart < 0)
    tlstart_idx = 0;

#pragma omp parallel default(none) shared(nvtxs, tlstart, tlstart_idx, tlsize, xadj, txadj, adjncy, tadjncy, adjwgt) private(nlocal, ej, ejend, ejstart, vi, ei, eiend) reduction(+:ntriangles) reduction(+:tnc)
  {
    int32_t *tmap1    = gk_i32smalloc(tlsize, -1, "tmap1");
    tmap1   -= tlstart; /* make indexing simpler */

    /* use tmap for the last tlsize rows */
    #pragma omp for schedule(dynamic, DYNAMIC_CHUNK)
    for (vj=tlstart_idx; vj<nvtxs; vj++) {
      if (1 || xadj[vj+1]-xadj[vj] < nvtxs-vj-1) {
        /* hash Adj(vj) */
        for (ej=ejstart=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++)
          tmap1[adjncy[ej]] = ej-ejstart;

        /* find intersections */
        for (ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
          vi = tadjncy[ej];
          for (nlocal=0, ei=xadj[vi]+tadjncy[ej+1], eiend=xadj[vi+1]; ei<eiend; ei++) {
            if (tmap1[adjncy[ei]] != -1) {
              #pragma omp atomic
              adjwgt[ei]++;
              #pragma omp atomic
              adjwgt[ejstart+tmap1[adjncy[ei]]]++;
              nlocal++;
            }
          }

          if (nlocal > 0) {
            ntriangles += nlocal;

            assert(adjncy[xadj[vi]+tadjncy[ej+1]-1] == vj);
            #pragma omp atomic
            adjwgt[xadj[vi]+tadjncy[ej+1]-1] += nlocal;
          }
        }

        /* reset tmap */
        for (ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++)
          tmap1[adjncy[ej]] = -1;
      }
      else { /* the row is dense */  /* TODO: This has not been updated */
        tnc++;
        /* find intersections */
        for (nlocal=0, ej=txadj[vj], ejend=txadj[vj+1]; ej<ejend; ej+=2) {
          vi = tadjncy[ej];
          nlocal += xadj[vi+1]-xadj[vi]-tadjncy[ej+1];
        }
        ntriangles += nlocal;
      }
    }

    tmap1 += tlstart;
    gk_free((void **)&tmap1, LTERM);
  }
  gk_stopwctimer(vault->timer_2);

  return ntriangles;
}
