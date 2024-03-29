
/******************************************************************************
 * INCLUDES
 *****************************************************************************/

#include "kt.h"
#include "kt_thread.h"
#include "kt_sbucket.h"
#include <stdbool.h>
#include <omp.h>
#include "timer.h"

/*
 * Size of chunk in dynamic for loops. Smaller values have been shown to be
 * very bad.
 */
#ifndef DYNAMIC_CHUNK
#define DYNAMIC_CHUNK 16
#endif

/*
 * Minimum size of frontier to bother with parallelizing the computation.
 */
#ifndef MIN_PAR_SIZE
#define MIN_PAR_SIZE 1
#endif


#define TIMER_PADDING 8
static double frontier_times[KT_MAX_THREADS * TIMER_PADDING] = { 0. };
static double update_times[KT_MAX_THREADS * TIMER_PADDING]   = { 0. };




/******************************************************************************
 * TYPES
 *****************************************************************************/


/*
 * Modified non-zero in CSR structures, but include inc/dec to act as a
 * doubly-linked list.
 */
typedef struct
{
  int32_t vj;
  int32_t dec;
  int32_t inc;
} aii_s;

/*
 * Modified rowptr from CSR structures, but includes start/end that we shrink.
 */
typedef struct
{
  int64_t start;
  int64_t end;
} xaii_s;

/*
 * Dynamic graph structure.
 */
typedef struct
{
  int32_t nvtxs;
  int64_t nedges; /* number of unique edges */

  int64_t total_support;
  int32_t * supports;

  /* maps aii[] space into edge_t[] space */
  ssize_t * ids;

  xaii_s * xaii;
  aii_s  * aii;
} dyn_graph_t;





/******************************************************************************
 * PRIVATE FUNCTIONS
 *****************************************************************************/

static dyn_graph_t * p_dgraph_alloc(
    int32_t const nvtxs,
    int64_t const nedges)
{
  dyn_graph_t * d_graph = gk_malloc(sizeof(*d_graph), "d_graph");

  d_graph->nvtxs  = nvtxs;
  d_graph->nedges = nedges;

  d_graph->xaii = gk_malloc((nvtxs+1) * sizeof(*d_graph->xaii), "xaii");
  d_graph->aii  = gk_malloc((2*nedges+1) * sizeof(*d_graph->aii), "aii");
  d_graph->ids  = gk_malloc((2*nedges+1) * sizeof(*d_graph->ids), "ids");
  d_graph->supports = gk_malloc(nedges * sizeof(*d_graph->supports), "supports");

  /* Go ahead and memset them in parallel to establish NUMA placement */
  par_memset(d_graph->xaii,  0, (nvtxs + 1) * sizeof(*d_graph->xaii));
  par_memset(d_graph->aii,   0, (2 * nedges + 1) * sizeof(*d_graph->aii));
  par_memset(d_graph->supports, 0, nedges * sizeof(*d_graph->supports));

  size_t bytes = 0;
  bytes += (nvtxs + 1) * sizeof(*d_graph->xaii);
  bytes += (2 * nedges + 1) * sizeof(*d_graph->aii);
  bytes += (2 * nedges + 1) * sizeof(*d_graph->ids);
  bytes += nedges * sizeof(*d_graph->supports);

  printf("DYN-GRAPH-BYTES: %0.2fGB\n",
      (double) bytes / (1024 * 1024 * 1024));

  return d_graph;
}

static void p_dgraph_free(
    dyn_graph_t * d_graph)
{
  gk_free((void **) &d_graph->aii, LTERM);
  gk_free((void **) &d_graph->xaii, LTERM);
  gk_free((void **) &d_graph->ids, LTERM);
  gk_free((void **) &d_graph->supports, LTERM);
  gk_free((void **) &d_graph, LTERM);
}




static void p_init_dgraph(
    dyn_graph_t * dgraph,
    gk_graph_t const * const ugraph,
    edge_t * const edges)
{
  /* extract data */
  int32_t const nvtxs = ugraph->nvtxs;
  ssize_t const * const restrict xadj = ugraph->xadj;
  int32_t const * const restrict adjncy = ugraph->adjncy;
  int32_t const * const restrict adjwgt = ugraph->iadjwgt;

  xaii_s * const restrict xaii = dgraph->xaii;
  aii_s  * const restrict aii = dgraph->aii;
  ssize_t * const restrict ids = dgraph->ids;
  int32_t * const restrict supports = dgraph->supports;

  /* initialize the start of each adj. list */
  #pragma omp parallel for
  for(int32_t vi=0; vi < nvtxs; ++vi) {
    dgraph->xaii[vi].start = 0;
  }

  /* Determine size of each adj. list */
  int64_t edge_ptr = 0;
  for(int32_t vi=0; vi < nvtxs; ++vi) {
    for(ssize_t ei = xadj[vi]; ei < xadj[vi+1]; ++ei) {
      if (adjwgt[ei] == 0) {
        continue;
      }

      xaii[vi].start++;
      xaii[adjncy[ei]].start++;

      edges[edge_ptr].vi = vi;
      edges[edge_ptr].vj = adjncy[ei];
      supports[edge_ptr] = adjwgt[ei];
      edge_ptr++;
    }
  }

  /* the MAKECSR equivalent */
  for(int32_t vi=1; vi < nvtxs; ++vi) {
    xaii[vi].start += xaii[vi-1].start;
  }
  for(int32_t vi=nvtxs; vi > 0; --vi) {
    xaii[vi].start = xaii[vi-1].start;
  }
  xaii[0].start = 0;

  /* populate dgraph two steps to ensure that sorted order is maintained */
  edge_ptr = 0;
  for(int32_t vi=0; vi < nvtxs; ++vi) {
    for(ssize_t ei = xadj[vi]; ei < xadj[vi+1]; ++ei) {
      if (adjwgt[ei] == 0) {
        continue;
      }

      int32_t const vj = adjncy[ei];
      aii[xaii[vj].start].vj  = vi;
      aii[xaii[vj].start].inc = 1;
      aii[xaii[vj].start].dec = 1;
      ids[xaii[vj].start]     = edge_ptr;

      edges[edge_ptr].eji     = xaii[vj].start++;
      edge_ptr++;
    }
  }
  edge_ptr = 0;
  for(int32_t vi=0; vi < nvtxs; ++vi) {
    for(ssize_t ei = xadj[vi]; ei < xadj[vi+1]; ++ei) {
      if (adjwgt[ei] == 0) {
        continue;
      }

      aii[xaii[vi].start].vj  = adjncy[ei];
      aii[xaii[vi].start].inc = 1;
      aii[xaii[vi].start].dec = 1;
      ids[xaii[vi].start]     = edge_ptr;

      edges[edge_ptr].eij     = xaii[vi].start++;
      edge_ptr++;
    }
  }
  /* the SHIFTCSR equivalent */
  for(int32_t vi=nvtxs; vi > 0; --vi) {
    xaii[vi].start = xaii[vi-1].start;
  }
  xaii[0].start = 0;

  /* record the end in xaii[vi] and from now own, you will be using that */
  for(int32_t vi=0; vi < nvtxs; ++vi) {
    xaii[vi].end = xaii[vi+1].start;
  }
}





static inline void p_intersect_lists(
    dyn_graph_t * const dgraph,
    int32_t const vi,
    int32_t const vj,
    thread_ws * ws,
    int32_t const * const supports,
    int64_t const edge_id,
    int32_t const min_support)
{
  assert(vi < vj);

  int32_t num_triangles = supports[edge_id];

  xaii_s const * const restrict xaii = dgraph->xaii;
  aii_s const * const restrict aii   = dgraph->aii;
  ssize_t const * const restrict ids = dgraph->ids;


  int64_t ei      = xaii[vi].end-1;
  int64_t ej      = xaii[vj].end-1;

  int64_t const eistart = xaii[vi].start;
  int64_t const ejstart = xaii[vj].start;

  /* decrease the support of the intersection */
  while (ei >= eistart && ej >= ejstart) {
    int32_t const vik = aii[ei].vj;
    int32_t const vjk = aii[ej].vj;
    if (vik > vjk) {
      ei  -= aii[ei].dec;
    } else if (vjk > vik) {
      ej  -= aii[ej].dec;
    } else {
      /* intersection! */

      /*
       * We need to queue the triangle for support updates. But, in order to
       * avoid duplicate triangle generation, we look at the other two edges
       * and ensure that we only queue the lowest-edge which is in the
       * frontier.
       *
       * Thus, we queue this triangle if it is the lowest-ID edge OR if the
       * lower-ID edges are not also in the frontier.
       */
      int32_t const edge_ei = ids[ei];
      int32_t const edge_ej = ids[ej];
      int64_t const min_edge         = gk_min(edge_ei, edge_ej);
      int32_t const min_edge_support = supports[min_edge];
      int64_t min_deleted_edge = edge_id;
      if((edge_ei < edge_id) && (supports[edge_ei] < min_support)) {
        min_deleted_edge = edge_ei;
      } else if((edge_ej < edge_id) && (supports[edge_ej] < min_support)) {
        min_deleted_edge = edge_ej;
      }
      /* if we should communicate */
      if(edge_id == min_deleted_edge) {
        triangle_t tri;

        /* sorting network so that u < v < w */
        tri.u = gk_min(vi, vik);
        tri.v = gk_min(gk_max(vi, vik), vj);
        tri.w = gk_max(vik, vj);

        tri.uv = edge_ei;
        tri.vw = edge_ej;

        int const dest = map_vtx_to_bucket(tri.u, ws);
        send_thread_tri_msg(&tri, dest, ws);
      }

      /* exit if we found all of the triangles */
      if(--num_triangles == 0) {
        break;
      }

#if VERBOSE
      printf("  TRIANGLE (%d %d %d)\n", 1+tri.u, 1+tri.v, 1+tri.w);
#endif

      ei  -= aii[ei].dec;
      ej  -= aii[ej].dec;
    }
  }
}


static inline void p_delete_edge(
    xaii_s * const restrict xaii,
    aii_s  * const restrict aii,
    int32_t const vi,
    int64_t const e_id)
{
  /* forward link */
  if(e_id == xaii[vi].start) {
    /* update start of list */
    xaii[vi].start += aii[e_id].inc;
  } else {
    /* delete node, re-route previous to next */
    aii[e_id - aii[e_id].dec].inc += aii[e_id].inc;
  }

  /* backwards link */
  if(e_id == xaii[vi].end - 1) {
    /* update end of list */
    xaii[vi].end -= aii[e_id].dec;
  } else {
    /* delete node, re-route previous to next */
    aii[e_id + aii[e_id].inc].dec += aii[e_id].dec;
  }
}



static void p_serial_peel(
    dyn_graph_t * dgraph,
    edge_t * const edges,
    int32_t const ktruss,
    int64_t * const restrict frontier_buf,
    int64_t const frontier_size,
    support_bucket_t * sbuckets,
    thread_ws * * thd_ws)
{

}




static int64_t p_gen_frontier(
    dyn_graph_t * dgraph,
    edge_t * const edges,
    int32_t const ktruss,
    int64_t * const restrict frontier_buf,
    int64_t const frontier_size,
    thread_ws * * thd_ws)
{
  xaii_s * const restrict xaii = dgraph->xaii;
  aii_s  * const restrict aii = dgraph->aii;
  int32_t * const restrict supports = dgraph->supports;

  int32_t const min_support = ktruss - 2;

  int64_t delta_support = 0;
  int64_t delta_edges = 0;

  #pragma omp parallel reduction(+: delta_edges, delta_support) \
      if(frontier_size >= MIN_PAR_SIZE)
  {
    int const tid = omp_get_thread_num();
    double my_timer;
    gk_clearwctimer(my_timer);
    gk_startwctimer(my_timer);

    /* Find edges to delete and discover triangles */
    #pragma omp for schedule(dynamic, DYNAMIC_CHUNK) nowait
    for(int64_t e_ptr = 0; e_ptr < frontier_size; ++e_ptr) {
      /* grab edge ID */
      int64_t const e = frontier_buf[e_ptr];
#if 0
      printf("  deleting %zd\n", e);
#endif
      assert(supports[e] < min_support);

      TimerBegin("p_intersect_lists", tid);
      if(supports[e] > 0) {
        int32_t const vi = edges[e].vi;
        int32_t const vj = edges[e].vj;
        assert(vi < vj);
        p_intersect_lists(dgraph, vi, vj, thd_ws[tid], supports, e,
            min_support);
      }
      TimerEnd("p_intersect_lists", tid);

      /* Send both sides of the edge. */
      TimerBegin("send_thread_epair_msg", tid);
      gk_epair_t msg;
      msg.key  = edges[e].vi;
      msg.val = edges[e].eij;
      send_thread_epair_msg(&msg,
          map_vtx_to_bucket(msg.key, thd_ws[tid]),
          thd_ws[tid]);

      msg.key = edges[e].vj;
      msg.val = edges[e].eji;
      send_thread_epair_msg(&msg,
          map_vtx_to_bucket(msg.key, thd_ws[tid]),
          thd_ws[tid]);

      ++delta_edges;
      delta_support += supports[e];

      /* (k-1) was the last valid k-truss for this edge */
      supports[e] = -(ktruss - 1);
      TimerEnd("send_thread_epair_msg", tid);
    } /* foreach edge */

    /* ensure all edges have been queued */
    gk_stopwctimer(my_timer);
    TimerBegin("barrier", tid);
    #pragma omp barrier
    gk_startwctimer(my_timer);
    TimerEnd("barrier", tid);

    /* actually delete edges from structure */
    int const num_buckets = thd_ws[tid]->num_buckets;
    #pragma omp for schedule(dynamic, 1) nowait
    for(int bucket=0; bucket < num_buckets; ++bucket) {
      int64_t num_del_edges = 0;
      TimerBegin("get_incoming_epair_bucket", tid);
      gk_epair_t * const del_pairs =
          get_incoming_epair_bucket(thd_ws, bucket, &num_del_edges);
      TimerEnd("get_incoming_epair_bucket", tid);

      /* sort edges to improve locality */
#if 0
      TimerBegin("gk_epairsorti", tid);
      gk_epairsorti(num_del_edges, del_pairs);
      TimerEnd("gk_epairsorti", tid);
#endif

      TimerBegin("p_delete_edge", tid);
      for(int32_t m=0; m < num_del_edges; ++m) {
        int64_t const e_id = del_pairs[m].val;
        int32_t const vi = del_pairs[m].key;
        p_delete_edge(xaii, aii, vi, e_id);

#if VERBOSE
        int32_t const vj = edges[e_id].vj;
        printf("  deleting %zd = %d (%d %d)\n",
            e_id, 1+aii[e_id].vj, 1 + vi, 1 + vj);
#endif
      }
      TimerEnd("p_delete_edge", tid);
    } /* foreach bucket */

    gk_stopwctimer(my_timer);
    frontier_times[omp_get_thread_num() * TIMER_PADDING] = my_timer;
  } /* omp parallel */


  dgraph->total_support -= delta_support;

  return delta_edges;
}



static int64_t p_update_supports(
    dyn_graph_t * dgraph,
    edge_t * edges,
    int32_t const ktruss,
    support_bucket_t * sbuckets,
    thread_ws * * thd_ws,
    int64_t const frontier_size)
{
  int const num_buckets = thd_ws[0]->num_buckets;

  int32_t * const restrict supports = dgraph->supports;

  int64_t nchanges = 0;
  #pragma omp parallel reduction(+: nchanges) \
      if(frontier_size > MIN_PAR_SIZE)
  {
    int const tid      = omp_get_thread_num();
    int const nthreads = omp_get_num_threads();

    double my_timer;
    gk_clearwctimer(my_timer);
    gk_startwctimer(my_timer);

    TimerBegin("update_1st", tid);
    /* go over incoming messages coming from each thread */
    #pragma omp for schedule(dynamic, 1) nowait
    for(int bucket=0; bucket < num_buckets; ++bucket) {
      /* grab bucket data */
      int32_t num_triangles = 0;
      TimerBegin("get_incoming_triangle_bucket", tid);
      triangle_t * const triangles =
          get_incoming_triangle_bucket(thd_ws, bucket, &num_triangles);
      TimerEnd("get_incoming_triangle_bucket", tid);

      TimerBegin("num_triangles_loop", tid);
      for(int t=0; t < num_triangles; ++t) {
        triangle_t const * const tri = &(triangles[t]);
#if VERBOSE
        printf("processing triangle (%d %d %d) = %zd %zd %zd\n",
            1 + tri->u, 1 + tri->v, 1 + tri->w, tri->uv, tri->vw);
#endif
        assert(tri->u < tri->v);
        assert(tri->v < tri->w);
        assert(tri->uv < dgraph->nedges);
        assert(tri->vw < dgraph->nedges);

        /* decrement first edge */
        int64_t const edge_id = tri->uv;
        if(supports[edge_id] > 0) {
          supports[edge_id]--;

          if(map_edge_to_bucket(edge_id, thd_ws[tid]) == bucket) {
            sbucket_update_edge(&(sbuckets[bucket]), edge_id, supports[edge_id],
                ktruss);
          } else {
            /* prepare to update support-based bucket */
            int const edge_bucket = map_edge_to_bucket(edge_id, thd_ws[tid]);
            send_thread_edge_msg(edge_id, edge_bucket, thd_ws[tid], 1);
          }
          ++nchanges;
#if VERBOSE
          printf("  (%zd) %d -> %d\n",
              tri->uv, supports[edge_id]+1, supports[edge_id]);
#endif
        }

        /* process second edge */
        int const bucket_dest = map_vtx_to_bucket(tri->w, thd_ws[tid]);
        if(bucket_dest == tri->w) {
          int64_t const edge_id = tri->vw;
          if(supports[edge_id] > 0) {
            supports[edge_id]--;

            if(map_edge_to_bucket(edge_id, thd_ws[tid]) == bucket) {
              sbucket_update_edge(&(sbuckets[bucket]), edge_id, supports[edge_id],
                  ktruss);
            } else {
              /* prepare to update support-based bucket */
              int const edge_bucket = map_edge_to_bucket(edge_id, thd_ws[tid]);
              send_thread_edge_msg(edge_id, edge_bucket, thd_ws[tid], 1);
            }
            ++nchanges;
          }

        } else {
          /* communicate and process in the next pass */
          send_thread_edge_msg(tri->vw, bucket_dest, thd_ws[tid], 0);
        }
      }
      TimerEnd("num_triangles_loop", tid);
    } /* foreach bucket */
    TimerEnd("update_1st", tid);

    gk_stopwctimer(my_timer);
    #pragma omp barrier
    gk_startwctimer(my_timer);


    TimerBegin("update_2nd", tid);
    #pragma omp for schedule(dynamic, 1) nowait
    for(int bucket=0; bucket < num_buckets; ++bucket) {
      int64_t num_edges = 0;
      TimerBegin("get_incoming_edge_bucket", tid);
      int64_t * const edges =
          get_incoming_edge_bucket(thd_ws, bucket, &num_edges, 0);
      TimerEnd("get_incoming_edge_bucket", tid);

      TimerBegin("e_loop", tid);
      for(int64_t e=0; e < num_edges; ++e) {
#if VERBOSE
        printf("thread %d received %zd edge-msgs from bucket %d\n",
            tid, num_edges, bucket);
#endif

        int64_t const edge_id = edges[e];
        if(supports[edge_id] > 0) {
          supports[edge_id]--;

          /* send to next stage */
          if(map_edge_to_bucket(edge_id, thd_ws[tid]) == bucket) {
            sbucket_update_edge(&(sbuckets[bucket]), edge_id, supports[edge_id],
                ktruss);
          } else {
            int const edge_bucket = map_edge_to_bucket(edge_id, thd_ws[tid]);
            send_thread_edge_msg(edge_id, edge_bucket, thd_ws[tid], 1);
          }
          ++nchanges;

#if VERBOSE
          printf("  (%zd) %d -> %d\n",
              edges[e], supports[edges[e]]+1, supports[edges[e]]);
#endif
        }
      }
      TimerEnd("e_loop", tid);
    } /* foreach bucket */
    TimerEnd("update_2nd", tid);


    gk_stopwctimer(my_timer);
    #pragma omp barrier
    gk_startwctimer(my_timer);

    /* finally, update all support-based buckets */
    TimerBegin("update_3rd", tid);
    #pragma omp for schedule(dynamic, 1) nowait
    for(int bucket=0; bucket < num_buckets; ++bucket) {
#if 1
      int64_t num_edges = 0;
      TimerBegin("get_incoming_edge_bucket2", tid);
      int64_t * const edges =
          get_incoming_edge_bucket(thd_ws, bucket, &num_edges, 1);
      TimerEnd("get_incoming_edge_bucket2", tid);

#if 0
      /* sort to improve locality */
      TimerBegin("u_gk_i64sorti", tid);
      gk_i64sorti(num_edges, edges);
      TimerEnd("u_gk_i64sorti", tid);
#endif

      /* foreach edge in bucket */
      for(int64_t e=0; e < num_edges; ++e) {
        int64_t const edge_id = edges[e];
        sbucket_update_edge(&(sbuckets[bucket]), edge_id, supports[edge_id],
            ktruss);
      }
#endif
    } /* foreach bucket */
    TimerEnd("update_3rd", tid);

    gk_stopwctimer(my_timer);
    update_times[omp_get_thread_num() * TIMER_PADDING] = my_timer;
  } /* end omp parallel */

#if VERBOSE
  printf("supports after update:");
  for(int64_t e=0; e < dgraph->nedges; ++e) {
    printf(" %d", supports[e]);
  }
  printf("\n");
#endif

  return nchanges;
}



/******************************************************************************
 * PUBLIC FUNCTIONS
 *****************************************************************************/

int64_t kt_twopass5b(params_t *params, vault_t *vault)
{
  InitTSC();
  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);

  /* Pull these out to save typing. */
  int32_t const nvtxs  = vault->ugraph->nvtxs;
  ssize_t const * const restrict xadj   = vault->ugraph->xadj;
  int32_t const * const restrict adjncy = vault->ugraph->adjncy;

  /* where the support values will be stored */
  vault->ugraph->iadjwgt = gk_i32malloc(xadj[nvtxs], "adjwgt");
  int32_t * const iadjwgt = vault->ugraph->iadjwgt;
  par_memset(iadjwgt, 0, xadj[nvtxs] * sizeof(*iadjwgt));
  gk_stopwctimer(vault->timer_tcsetup);

  gk_startwctimer(vault->timer_esupport);
  int64_t ntriangles = kt_ComputeEdgeSupportPar(params, vault);
  gk_stopwctimer(vault->timer_esupport);
  gk_graph_Free(&(vault->lgraph));

  gk_startwctimer(vault->timer_ktsetup);
  /*
   * Allocate/initialize edges. We only bother with edges that have non-zero
   * support.
   */
  int64_t nedges = count_nnz(xadj[nvtxs], iadjwgt);
  edge_t * edges = gk_malloc((nedges+1) * sizeof(*edges), "edges");
  par_memset(edges, 0, (nedges + 1) * sizeof(*edges));

  /* Allocate and fill the dynamic graph */
  dyn_graph_t * dgraph = p_dgraph_alloc(nvtxs, nedges);
  p_init_dgraph(dgraph, vault->ugraph, edges);

  gk_stopwctimer(vault->timer_ktsetup);

  size_t edge_bytes = (nedges+1) * sizeof(*edges);
  printf("EDGES-BYTES: %0.3fGB\n",
      (double) edge_bytes / (1024. * 1024. * 1024));


  printf("#triangles before peeling: %zd\n", ntriangles);
  ntriangles = 0;


  int64_t edges_left = nedges;
  dgraph->total_support = count_support(nedges, dgraph->supports);

  printf("THREADS: %d\n", omp_get_max_threads());
  printf("DYNAMIC_CHUNK: %d\n", DYNAMIC_CHUNK);
  printf("KT_BUCKETS_PER_THREAD: %d\n", KT_BUCKETS_PER_THREAD);
  printf("MIN_PAR_SIZE: %d\n", MIN_PAR_SIZE);

  printf("\n");

  thread_ws * * thd_ws = alloc_thread_ws_big(vault->graph);

  /* Setup support-based buckets. */
  support_bucket_t * support_buckets = sbucket_alloc(edges, dgraph->supports,
      dgraph->nedges, thd_ws);

  /* XXX should be much smaller and resize.. */
  int64_t * frontier_buf = gk_malloc(dgraph->nedges * sizeof(*frontier_buf),
      "frontier_buf");
  /*
   * Main loop.
   */
  /* initial status */
  printf("k: %7d; edges-left: %7zd (%6.2f%%), total-support: %7zd, time (s): %6.3f\n",
      3, edges_left, 100. * (double)edges_left / (double)nedges,
      dgraph->total_support,
      0.);

  /* timers */
  double timer_currk;
  double timer_frontier;
  double timer_updates;
  double total_frontier_time;
  double total_update_time;
  gk_clearwctimer(timer_currk);
  gk_clearwctimer(timer_frontier);
  gk_clearwctimer(timer_updates);
  gk_clearwctimer(total_frontier_time);
  gk_clearwctimer(total_update_time);


  gk_startwctimer(timer_currk);
  gk_startwctimer(vault->timer_ktpeeling);
  int32_t ktruss = 4;
  while(edges_left > 0) {

    gk_startwctimer(timer_frontier);
    /* ktruss-3 is everything we need to remove */
    int64_t const frontier_size = sbucket_get_frontier(support_buckets,
        ktruss-3, frontier_buf);
    int64_t const delta_edges = p_gen_frontier(dgraph, edges, ktruss,
        frontier_buf, frontier_size, thd_ws);
    gk_stopwctimer(timer_frontier);

    int64_t delta_support;
    gk_startwctimer(timer_updates);
    if(delta_edges < edges_left) {
      delta_support = p_update_supports(dgraph, edges, ktruss,
          support_buckets, thd_ws, frontier_size);
    } else {
      /* if we have removed the remaining edges */
      delta_support = dgraph->total_support;
      dgraph->total_support = 0;
    }
    gk_stopwctimer(timer_updates);

#if 0
    printf("    frontier:");
    thread_time_stats(frontier_times, omp_get_max_threads(), TIMER_PADDING);
    printf("    supports:");
    thread_time_stats(update_times, omp_get_max_threads(), TIMER_PADDING);
#endif

    edges_left -= delta_edges;
    dgraph->total_support -= delta_support;

    if(edges_left == 0 || (delta_edges == 0 && delta_support == 0)) {
      gk_stopwctimer(timer_currk);
      printf("k: %7d; edges-left: %7zd (%6.2f%%), total-support: %7zd, time (s): %6.3f\n",
          ktruss, edges_left, 100. * (double)edges_left / (double)nedges,
          dgraph->total_support,
          timer_currk);
      printf("  frontier: %6.3fs, updates: %6.3fs\n",
          timer_frontier, timer_updates);

      total_frontier_time += timer_frontier;
      total_update_time += timer_updates;

      ++ktruss;
      gk_clearwctimer(timer_frontier);
      gk_clearwctimer(timer_updates);
      gk_clearwctimer(timer_currk);

      gk_startwctimer(timer_currk);
    }
  } /* end main loop */
  gk_stopwctimer(vault->timer_ktpeeling);

  printf("\n");
  printf("FRONTIER: %0.3fs UPDATE: %0.3fs\n",
      total_frontier_time, total_update_time);


  printf("#triangles after peeling: %zd\n", ntriangles);

  TimerPrint();

  /* adjust for George's counting :-) */
  if(params->outfile != NULL) {
    #pragma omp parallel for schedule(static)
    for(int64_t e=0; e < nedges; ++e) {
      dgraph->supports[e] += 2;
    }

    /* create the output of the decomposition */
    kt_Sups2KTEdges(params, vault, ktruss-1, dgraph->supports);
  }

  /*
   * Cleanup.
   */
  p_dgraph_free(dgraph);
  gk_free((void **)&edges, LTERM);
  gk_free((void **)&frontier_buf, LTERM);
  free_thread_ws(thd_ws);
  sbucket_free(support_buckets);

  return ntriangles;
}


