/*!
\file
\brief Data structures used in the program
\date Started 5/10/2017
\author George
\version\verbatim $Id: struct.h 21155 2017-06-09 21:36:55Z karypis $ \endverbatim
*/

#ifndef _STRUCT_KT_H_
#define _STRUCT_KT_H_


typedef struct
{
  int32_t vi; /* source vtx */
  int32_t vj; /* dest vtx */

  ssize_t eij; /* index of edge (vi, vj) */
  ssize_t eji; /* index of edge (vj, vi) */
} edge_t;


/*************************************************************************
* an edge in the k-truss decomposition
**************************************************************************/
typedef struct {
  int32_t vi, vj, k;
} ktedge_t;

/*************************************************************************
* the data vault
**************************************************************************/
typedef struct {
  gk_graph_t *graph;      /* the graph */
  gk_graph_t *ugraph;     /* the upper triangular graph */
  gk_graph_t *lgraph;     /* the lower triangular graph (with potentially indexing
                             information in U */

  int32_t *perm;        /* perm[old-vtx-num]  => new-vtx-num */
  int32_t *iperm;       /* iperm[new-vtx-num] => old-vtx-num */

  int64_t nedges;       /* the number of edges in ktedges */
  ktedge_t *ktedges;    /* the k-truss decomposition */
  int32_t ktmax;        /* the maximum k found in k-truss */

  /* timers */
  double timer_global;
  double timer_io;
  double timer_tcsetup;
  double timer_ktsetup;
  double timer_esupport;
  double timer_ktpeeling;
  double timer_1;
  double timer_2;
  double timer_3;
  double timer_4;
  double timer_5;
  double timer_6;

  /* Opttional parameters for streaming scenario */
  int32_t *kt;
  ssize_t *endIdx;
  ssize_t *revIdx;

} vault_t;


/*************************************************************************
* run parameters
**************************************************************************/
typedef struct {
  int kttype;           /* The algorithm to use for k-truss */
  int iftype;           /* The format of the input file */

  int jbsize;           /* The size of the J block in words (query block) */
  int ibsize;           /* The size of the I block in words (database block) */

  int seed;             /* Seed for the random number generator */
  int dbglvl;           /* Debuging information */

  char *infile;         /* The file storing the input data */
  char *outfile;        /* The file storing the output data */

  int stream;           /* Boolean value to denote whether accept new edges */
} params_t;



/******************************************************************************
 * SHADEN TYPES
 *****************************************************************************/
typedef struct
{
  /* u < v < w */
  int32_t u;
  int32_t v;
  int32_t w;
  /* edge indices */
  int64_t uv;
  //int64_t uw;
  int64_t vw;
} triangle_t;


typedef struct
{
  int32_t v;
  int32_t w;
  int64_t vw;
} pair_t;


/*************************************************************************
* Data structures for streaming k-truss
**************************************************************************/

typedef struct eedge_t eedge_t;
struct eedge_t {
  int32_t v; /* 4 bytes */
  eedge_t *revPtr; /* 8 bytes */
  int8_t k, td, visited; /* 1*3 bytes */
  int8_t mtd; 
}; /* 4 + 8 + 1*5 = 17 bytes */

typedef struct {
  int32_t nvtxs;
  ssize_t *xadj;
  eedge_t *edges;
} ggraph_t;

typedef struct edgeblock_t edgeblock_t;
struct edgeblock_t{
  int32_t src; /* 4 bytes */
  /* Change below to an array to include multiple edges in a block? */
  eedge_t edge; /* 13 bytes */
  edgeblock_t *next; /* 8 bytes */
}; /* 25 bytes */

#endif
