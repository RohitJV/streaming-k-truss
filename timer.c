#include <float.h>
#include <omp.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "timer.h"

static int initialized = 0;
static int fifo_head = 0;
static int fifo[FIFO_DEPTH];
static int n_code_blocks = 0;
static uint64_t timer_in[MAX_THREADS][MAX_CODE_BLOCKS];
static double timer_sum[MAX_THREADS][MAX_CODE_BLOCKS];
static char code_block_label[MAX_CODE_BLOCKS][MAX_LABEL + 1];
static int block_starts[MAX_THREADS][MAX_CODE_BLOCKS];
static int block_stops[MAX_THREADS][MAX_CODE_BLOCKS];


static double g_ticks_persecond = 0.0;

void InitTSC(void) {
  uint64_t start_tick = ReadTSC();
  sleep(1);
  uint64_t end_tick = ReadTSC();

  g_ticks_persecond = (double)(end_tick - start_tick);

  memset(block_starts, 0, sizeof(int) * MAX_CODE_BLOCKS * MAX_THREADS);
  memset(block_stops, 0, sizeof(int) * MAX_CODE_BLOCKS * MAX_THREADS);
  memset(timer_sum, 0, sizeof(double) * MAX_CODE_BLOCKS * MAX_THREADS);

  initialized = 1;
}



double ElapsedTime(uint64_t ticks) {
  if (g_ticks_persecond == 0.0) {
    fprintf(stderr, "TSC timer has not been initialized.\n");
    return 0.0;
  } else {
    return (ticks / g_ticks_persecond);
  }
}

static int IndexFromLabel(char *this_label) {
  // check FIFO
  if (n_code_blocks > FIFO_DEPTH) {
    for (int k = 0; k < FIFO_DEPTH; k++) {
      int idx = fifo[k];
      if (0 == strncmp(code_block_label[idx], this_label, MAX_LABEL)) {
        return idx;
      }
    }
  }

  // check from the end
  for (int idx = n_code_blocks - 1; idx >= 0; idx--) {
    if (0 == strncmp(code_block_label[idx], this_label, MAX_LABEL)) {
      return idx;
    }
  }

  // not exist, a new code block
  strncpy(code_block_label[n_code_blocks], this_label, MAX_LABEL);
  code_block_label[n_code_blocks][MAX_LABEL] = '\0';
  fifo[fifo_head] = n_code_blocks;
  fifo_head = (fifo_head + 1) % FIFO_DEPTH;
  n_code_blocks++;
  return (n_code_blocks - 1);
}

/*
 * These are defined in timer.h if timers are disabled.
 */
#ifdef USE_TIMER_
void TimerBegin(char *this_label, int tid) {
  if (n_code_blocks == MAX_CODE_BLOCKS) {
    return;
  }

  assert(initialized == 1);

  // find this label
  int idx = IndexFromLabel(this_label);
  if (idx >= MAX_CODE_BLOCKS) {
    return;
  }

  block_starts[tid][idx] += 1;
  timer_in[tid][idx] = ReadTSC();
}

void TimerEnd(char *this_label, int tid) {
  uint64_t t_end = ReadTSC();
  if (n_code_blocks == MAX_CODE_BLOCKS) {
    return;
  }

  // find this label
  int idx = IndexFromLabel(this_label);
  block_stops[tid][idx] += 1;
  timer_sum[tid][idx] += ElapsedTime(t_end - timer_in[tid][idx]);
}

void TimerReset(void) {
  if (n_code_blocks > 0) {
    n_code_blocks = 0;
    fifo_head = 0;
    initialized = 0;
  }
}
#endif




void TimerPrint(void) {
#ifdef USE_TIMER_
  printf("\n");
  printf("---------------------------------------------------------------------"
         "---------------------------------------------------------------------"
         "------------\n");
  printf(
      "Perf.  summary:                  avg(secs)    min(secs)    max(secs)    imbalance(%%)\n");
  printf("---------------------------------------------------------------------"
         "---------------------------------------------------------------------"
         "------------\n");
  int nthreads = omp_get_max_threads();
  for (int j = 0; j < n_code_blocks; j++) {
    double min_val = DBL_MAX;
    double max_val = 0.0;
    double avg_val = 0.0;
    double imb_val = 0.0; /* imbalance */
    for (int tid = 0; tid < nthreads; tid++) {
      if (block_starts[tid][j] == block_stops[tid][j]) {
        // time
        double my_val = timer_sum[tid][j];
        min_val = (min_val > my_val ? my_val : min_val);
        max_val = (max_val < my_val ? my_val : max_val);
        avg_val += my_val;
      } else {
        printf("mismatch in starts/stops for code block '%s'\n",
               code_block_label[j]);
        printf("  starts = %d\n", block_starts[tid][j]);
        printf("  stops  = %d\n", block_stops[tid][j]);
        break;
      } // if (block_starts[j] == block_stops[j])
    }
    avg_val = avg_val / nthreads;
    imb_val = 100. * ((max_val - avg_val) / avg_val);
    printf("%-28s  %12.3f %12.3f %12.3f %12.3f%%\n", code_block_label[j], avg_val,
           min_val, max_val, imb_val);
  } // for (int j = 0; j< n_code_blocks; j++)
#endif
}
