/*
 *                             The MIT License
 *
 * This file is part of QuickEdit library.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef ALIGN_BENCHMARK_PARAMS_H_
#define ALIGN_BENCHMARK_PARAMS_H_

#include "utils/commons.h"
#include "benchmark/benchmark_utils.h"

// Enable external benchmarks
#define EXTERNAL_BENCHMARKS

/*
 * Algorithms
 */
typedef enum {
  // Test
  alignment_test,
  // Edit
  alignment_edit_dp,
  alignment_edit_dp_banded,
  //
  alignment_edit_bpm,
  alignment_edit_bpm_banded,
  alignment_edit_bpm_banded_unaligned,
  alignment_edit_bpm_banded_blocking,
  alignment_edit_bpm_windowed,
  alignment_edit_bpm_quicked,
#ifdef EXTERNAL_BENCHMARKS
  // External algorithms
  alignment_bitpal_edit,
  alignment_daligner,
  alignment_diffutils,
  alignment_edlib,
  alignment_lv89,
  alignment_parasail_nw_stripped,
  alignment_parasail_nw_scan,
  alignment_parasail_nw_diag,
  alignment_parasail_nw_banded,
  alignment_scrooge,
  alignment_seqan_edit,
  alignment_seqan_edit_bpm,
#endif
} alignment_algorithm_type;

/*
 * Align-benchmark Parameters
 */
typedef struct {
  // Algorithm
  alignment_algorithm_type algorithm;
  // I/O
  char *input_filename;
  char *output_filename;
  bool output_full;
  // I/O internals
  FILE* input_file;
  char* line1;
  char* line2;
  size_t line1_allocated;
  size_t line2_allocated;
  FILE* output_file;
  // Alignment form
  bool endsfree;
  double pattern_begin_free;
  double text_begin_free;
  double pattern_end_free;
  double text_end_free;
  // Other algorithms parameters
  int bandwidth;
  int window_size;
  int overlap_size;
  window_config_t window_config;
#ifdef EXTERNAL_BENCHMARKS
  /* ... */
#endif
  // Misc
  bool check_display;
  bool check_correct;
  bool check_score;
  bool check_alignments;
  int check_bandwidth;
  int plot;
  // Profile
  profiler_timer_t timer_global;
  // System
  int num_threads;
  int batch_size;
  int progress;
  int verbose;
} align_bench_params_t;

// Defaults
extern align_bench_params_t parameters;

/*
 * Menu
 */
void usage();

/*
 * Parse arguments
 */
void parse_arguments(
    int argc,
    char** argv);

#endif /* ALIGN_BENCHMARK_PARAMS_H_ */
