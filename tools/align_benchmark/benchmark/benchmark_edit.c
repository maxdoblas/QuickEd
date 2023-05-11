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

#include "../../../alignment/bpm.h"
#include "benchmark/benchmark_utils.h"
#include "benchmark/benchmark_check.h"
#include "edit/edit_dp.h"
#include "../../../alignment/bpm_windowed.h"

/*
 * Benchmark Edit
 */
void benchmark_edit_bpm(
    align_input_t* const align_input) {
  // Allocate
  bpm_pattern_t bpm_pattern;
  bpm_pattern_compile(
      &bpm_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  bpm_matrix_t bpm_matrix;
  bpm_matrix_allocate(
      &bpm_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  bpm_compute(
      &bpm_matrix,&bpm_pattern,align_input->text,
      align_input->text_length,align_input->pattern_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,bpm_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,bpm_matrix.cigar);
  }
  // Free
  bpm_pattern_free(&bpm_pattern,align_input->mm_allocator);
  bpm_matrix_free(&bpm_matrix,align_input->mm_allocator);
}
void benchmark_edit_bpm_banded(
    align_input_t* const align_input,
    const int bandwidth) {
  // Allocate
  bpm_pattern_t bpm_pattern;
  bpm_pattern_compile(
      &bpm_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  bpm_matrix_t bpm_matrix;
  bpm_matrix_allocate(
      &bpm_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  bpm_compute_banded(
      &bpm_matrix,&bpm_pattern,align_input->text,
      align_input->text_length,bandwidth,align_input->pattern_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,bpm_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,bpm_matrix.cigar);
  }
  // Free
  bpm_pattern_free(&bpm_pattern,align_input->mm_allocator);
  bpm_matrix_free(&bpm_matrix,align_input->mm_allocator);
}
void benchmark_edit_bpm_windowed(
    align_input_t* const align_input) {
  // Allocate
  windowed_pattern_t windowed_pattern;
  windowed_pattern_compile(
      &windowed_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  windowed_matrix_t windowed_matrix;
  windowed_matrix_allocate(
      &windowed_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  windowed_compute(
      &windowed_matrix,&windowed_pattern,align_input->text,
      align_input->text_length,align_input->pattern_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,windowed_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,windowed_matrix.cigar);
  }
  // Free
  windowed_pattern_free(&windowed_pattern,align_input->mm_allocator);
  windowed_matrix_free(&windowed_matrix,align_input->mm_allocator);
}
void benchmark_edit_dp(
    align_input_t* const align_input) {
  // Parameters
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  // Allocate
  score_matrix_t score_matrix;
  score_matrix_allocate(
      &score_matrix,pattern_length+1,
      text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      pattern_length+text_length);
  // Align
  timer_start(&align_input->timer);
  edit_dp_align(&score_matrix,
      align_input->pattern,pattern_length,
      align_input->text,text_length,cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,cigar);
  }
  // Free
  score_matrix_free(&score_matrix);
  cigar_free(cigar);
}
void benchmark_edit_dp_banded(
    align_input_t* const align_input,
    const int bandwidth) {
  // Parameters
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  // Allocate
  score_matrix_t score_matrix;
  score_matrix_allocate(
      &score_matrix,pattern_length+1,
      text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      pattern_length+text_length);
  // Align
  timer_start(&align_input->timer);
  edit_dp_align_banded(&score_matrix,
      align_input->pattern,pattern_length,
      align_input->text,text_length,
      bandwidth,cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,cigar);
  }
  // Free
  score_matrix_free(&score_matrix);
  cigar_free(cigar);
}
