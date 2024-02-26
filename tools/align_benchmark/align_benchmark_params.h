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

#include "utils/include/commons.h"
#include "benchmark/benchmark_utils.h"

/*
 * Algorithms
 */
typedef enum {
  // Edit
  alignment_edit_dp,
  alignment_edit_dp_banded,
  alignment_edit_bpm,
  // QuickEd lib
  alignment_edit_banded,
  alignment_edit_banded_hirschberg,
  alignment_edit_windowed,
  alignment_edit_quicked,
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
  // Other algorithms parameters
  int bandwidth;
  int window_size;
  int overlap_size;
  int hew_threshold;
  int hew_percentage;
  bool force_scalar;
  bool only_score;
  // Misc
  bool check_display;
  bool check_correct;
  bool check_score;
  bool check_alignments;
  int check_bandwidth;
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
void usage(void);

/*
 * Parse arguments
 */
void parse_arguments(
    int argc,
    char** argv);

#endif /* ALIGN_BENCHMARK_PARAMS_H_ */
