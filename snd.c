
#include "kt.h"
#include <stdbool.h>

#ifndef MAX_NTHREADS
#define MAX_NTHREADS 272
#endif


/******************************************************************************
 * PRIVATE FUNCTIONS
 *****************************************************************************/


/**
* @brief Map an edge id to the owning vertex (i.e., u in a (u,v) edge).
*
* @param adj The xadj list (i.e., row_ptr in a CSR).
* @param edge The edge index.
*
* @return Which vertex the edge belongs to.
*/
static int32_t p_lookup_vtx(
    ssize_t const * const adj,
    int32_t const edge)
{
  int32_t v = 0;
  while(adj[v+1] <= edge) {
    ++v;
  }

  assert(adj[v] <= edge);
  assert(adj[v+1] > edge);

  return v;
}


/**
* @brief Intersect two adjacency lists (`adj_u` and `adj_v`)
*
* @param adj_u The neighbors of vertex u.
* @param len_u The length of `adj_u`.
* @param adj_u_offset This is how far into the xadj we already are. This allows
*                     us to translate len_u into a global graph edge.
* @param adj_v The neighbors of vertex v.
* @param len_v The length of `adj_v`.
* @param[out] triangles An array of the edge IDs which complete the discovered
*                       triangles.
* @param max_triangles The maximum number of triangles which we should find in
*                      the intersection.
*
* @return The number of discovered triangles.
*/
static int32_t p_intersect_lists(
    int32_t * const restrict adj_u,
    ssize_t const len_u,
    ssize_t const adj_u_offset,
    int32_t * const restrict adj_v,
    ssize_t const len_v,
    int32_t * const restrict triangles,
    int32_t const max_triangles)
{
  if(max_triangles == 0) {
    return 0;
  }

  int32_t num_found = 0;

  /* Linear merge to find intersections. We go in reverse because high-degree
   * vertices are placed at the end, and are thus more likely to be found in
   * the intersections. */
  int32_t u_ptr = len_u - 1;
  int32_t v_ptr = len_v - 1;
  while((u_ptr >= 0) && (v_ptr >= 0)) {
    int32_t const u = adj_u[u_ptr];
    int32_t const v = adj_v[v_ptr];
    if(u < v) {
      --v_ptr;
    } else if(v < u) {
      --u_ptr;
    } else {
      triangles[num_found++] = u_ptr + adj_u_offset;
      if(num_found == max_triangles) {
        return num_found;
      }
      --u_ptr;
      --v_ptr;
    }
  }

  return num_found;
}


static void p_find_triangles(
    gk_graph_t const * const lgraph,
    gk_graph_t const * const ugraph,
    int32_t const num_triangles,
    int32_t const * const restrict supports,
    int32_t       * const restrict h_index,
    int32_t const u,
    int32_t const v)
{
  int32_t found_triangles = 0;

  /*
   * For each triangle, we need to find vertex 'W' which completes the
   * triangle. There are three cases to consider:
   *   (1) (u, v, W) -> 'W' will be in ugraph[u] and ugraph[v].
   *   (2) (u, W, v) -> 'W' will be in ugraph[u] and lgraph[v].
   *   (3) (W, u, v) -> 'W' will be in lgraph[u] and lgraph[v].
   */

  /* XXX add software prefetching of adj[u] and adj[w]? */

  int32_t nnbrs_u = ugraph->xadj[u+1] - ugraph->xadj[u];
  int32_t nnbrs_v = ugraph->xadj[v+1] - ugraph->xadj[v];
  int32_t * adj_u = &(ugraph->adjncy[ugraph->xadj[u]]);
  int32_t * adj_v = &(ugraph->adjncy[ugraph->xadj[v]]);

  /* (u, v, W) */
  if(found_triangles != num_triangles) {
    int32_t const new_triangles = p_intersect_lists(
        adj_u, nnbrs_u,
        ugraph->xadj[u],
        adj_v, nnbrs_v,
        &(h_index[found_triangles]), num_triangles - found_triangles);

    found_triangles += new_triangles;
  }

  /* (u, W, v) */
  if(found_triangles != num_triangles) {
    nnbrs_v = lgraph->xadj[v+1] - lgraph->xadj[v];
    adj_v   = &(lgraph->adjncy[lgraph->xadj[v]]);
    int32_t const new_triangles = p_intersect_lists(
        adj_u, nnbrs_u,
        ugraph->xadj[u],
        adj_v, nnbrs_v,
        &(h_index[found_triangles]), num_triangles - found_triangles);

    found_triangles += new_triangles;
  }


  /* (W, u, v) */
  if(found_triangles != num_triangles) {
    nnbrs_u = lgraph->xadj[u+1] - lgraph->xadj[u];
    adj_u   = &(lgraph->adjncy[lgraph->xadj[u]]);
    int32_t const new_triangles = p_intersect_lists(
        adj_u, nnbrs_u,
        lgraph->xadj[u],
        adj_v, nnbrs_v,
        &(h_index[found_triangles]), num_triangles - found_triangles);

    /* we have to translate the edges in lgraph(u) to ugraph(W) */
    for(int32_t t=0; t < new_triangles; ++t) {
      int32_t const source_vtx = lgraph->adjncy[h_index[found_triangles + t]];
      assert(source_vtx < u);

      /* now find u in xadj */
      ssize_t const start = ugraph->xadj[source_vtx];
      ssize_t const stop = ugraph->xadj[source_vtx+1];
      for(ssize_t e=stop; e >= start; --e) {
        if(ugraph->adjncy[e] == u) {
          h_index[found_triangles + t] = e;
        }
      }
    }
    found_triangles += new_triangles;
  }
  assert(found_triangles == num_triangles);

  /* translate edge IDs to actual support values */
  for(int32_t t=0; t < num_triangles; ++t) {
    assert(h_index[t] < ugraph->xadj[ugraph->nvtxs]);
    h_index[t] = supports[h_index[t]];
  }
}



static int32_t p_compute_hindex(
    int32_t const * const restrict vals,
    int32_t       * const restrict buffer,
    int32_t const N)
{
  for(int32_t i=0; i < N+1; ++i) {
    buffer[i] = 0;
  }

  for(int32_t i=0; i < N; ++i) {
    int32_t idx = 0;
    if(vals[i] < N) {
      idx = vals[i];
    } else {
      idx = N;
    }
    ++buffer[idx];
  }

  int32_t sum = 0;
  for(int32_t i=N; i >= 0; --i) {
    sum += buffer[i];
    if(sum >= i) {
      return i;
    }
  }

	assert(false);
  return -1;
}


static int32_t p_update_edge(
    gk_graph_t const * const lgraph,
    gk_graph_t const * const ugraph,
    int32_t const * const restrict supports,
    int32_t       * const restrict h_index,
    int32_t       * const restrict h_index_buf,
    int64_t const edge_idx)
{
  int32_t const u = p_lookup_vtx(ugraph->xadj, edge_idx);
  int32_t const v = ugraph->adjncy[edge_idx];
  int32_t const num_triangles = ugraph->iadjwgt[edge_idx];

  p_find_triangles(lgraph, ugraph, num_triangles, supports, h_index, u, v);

  return p_compute_hindex(h_index, h_index_buf, num_triangles);
}






/******************************************************************************
 * PUBLIC FUNCTIONS
 *****************************************************************************/

int64_t kt_snd(params_t *params, vault_t *vault)
{
  printf("THREADS: %d\n", omp_get_max_threads());

  /*
   * Grab upper and lower trangular portion of graph.
   */
  gk_startwctimer(vault->timer_tcsetup);
  vault->ugraph = kt_PreprocessAndExtractUpper(params, vault);
  vault->lgraph = kt_TransposeUforJIK(params, vault->ugraph);
  gk_stopwctimer(vault->timer_tcsetup);

  int32_t const nvtxs  = vault->ugraph->nvtxs;
  int64_t const nedges = vault->ugraph->xadj[nvtxs];

  /*
   * Compute initial supports and count the number of edges with support
   * greater than zero.
   */
  gk_startwctimer(vault->timer_esupport);
  vault->ugraph->iadjwgt = gk_i32malloc(nedges, "iadjwgt");
  par_memset(vault->ugraph->iadjwgt, 0, nedges * sizeof(*vault->ugraph->iadjwgt));
  int32_t * supports     = gk_i32malloc(nedges, "supports");
  int32_t * new_supports = gk_i32malloc(nedges, "new_supports");

  par_memset(supports, 0, nedges * sizeof(*supports));

  int64_t const ntriangles = kt_ComputeEdgeSupport(params, vault);
  par_memcpy(supports, vault->ugraph->iadjwgt, nedges * sizeof(*supports));

  int64_t const nz_edges = count_nnz(nedges, supports);
  gk_stopwctimer(vault->timer_esupport);

  printf("Found |V|=%d |E|=%ld |T|=%ld NZ-SUPPORTS=%ld (%0.1f%%)\n",
      nvtxs, nedges, ntriangles, nz_edges,
      100. * (double) nz_edges / (double) nedges);

  gk_graph_Free(&vault->lgraph);
  vault->lgraph = transpose_graph(vault->ugraph);

  int32_t const max_support = max_elem(supports, nedges);

  int32_t * h_index[MAX_NTHREADS];
  int32_t * h_index_buf[MAX_NTHREADS];
  #pragma omp parallel
  {
    int const tid = omp_get_thread_num();
    h_index[tid] = gk_malloc(max_support * sizeof(*h_index), "h_index");
    h_index_buf[tid] = gk_malloc((max_support+1) * sizeof(*h_index),
        "h_index_buf");
  }

  bool done = false;
  /*
   * Main loop.
   */
  gk_startwctimer(vault->timer_ktpeeling);
  while(!done) {
    int64_t nchanges = 0;
    #pragma omp parallel for schedule(dynamic, 256) reduction(+: nchanges)
    for(int64_t e=0; e < nedges; ++e) {
      int const tid = omp_get_thread_num();
      int32_t old_support = supports[e];
      supports[e] = p_update_edge(vault->lgraph, vault->ugraph,
          supports, h_index[tid], h_index_buf[tid], e);

      if(supports[e] != old_support) {
        ++nchanges;
      }
    }

    done = (nchanges == 0);

    printf("nchanges: %zd\n", nchanges);
  } /* end main loop */
  gk_stopwctimer(vault->timer_ktpeeling);


  /* cleanup thread data */
  #pragma omp parallel
  {
    int const tid = omp_get_thread_num();
    gk_free((void **) &h_index[tid], LTERM);
    gk_free((void **) &h_index_buf[tid], LTERM);
  }


  int32_t max_ktruss = 0;
  for(int64_t e=0; e < nedges; ++e) {
    /* +2 because of the k-truss definition... */
    supports[e] += 2;
    max_ktruss = gk_max(max_ktruss, supports[e]);
  }
  /* +2 because of the k-truss definition... */
  printf("\nMAX K-TRUSS: %d\n\n", max_ktruss);

  gk_free((void **) &new_supports, LTERM);

  return (int64_t) ntriangles;
}

