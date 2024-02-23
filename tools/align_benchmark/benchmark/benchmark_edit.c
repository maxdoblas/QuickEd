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

#include "benchmark/benchmark_utils.h"
#include "benchmark/benchmark_check.h"
#include "edit/edit_dp.h"
#include "../../../quicked/quicked.h"
#include "utils/include/commons.h"

#define BPM_W64_LENGTH UINT64_LENGTH


/*
 * Benchmark Edit
 */
//void benchmark_edit_bpm(
//    align_input_t* const align_input) {
//  // Allocate
//  bpm_pattern_t bpm_pattern;
//  bpm_pattern_compile(
//      &bpm_pattern,align_input->pattern,
//      align_input->pattern_length,align_input->mm_allocator);
//  bpm_matrix_t bpm_matrix;
//  bpm_matrix_allocate(
//      &bpm_matrix,align_input->pattern_length,
//      align_input->text_length,align_input->mm_allocator);
//  // Align
//  timer_start(&align_input->timer);
//  bpm_compute(
//      &bpm_matrix,&bpm_pattern,align_input->text,
//      align_input->text_length,align_input->pattern_length);
//  timer_stop(&align_input->timer);
//  // DEBUG
//  if (align_input->debug_flags) {
//    benchmark_check_alignment(align_input,bpm_matrix.cigar);
//  }
//  // Output
//  if (align_input->output_file) {
//    benchmark_print_output(align_input,false,bpm_matrix.cigar);
//  }
//  // Free
//  bpm_pattern_free(&bpm_pattern,align_input->mm_allocator);
//  bpm_matrix_free(&bpm_matrix,align_input->mm_allocator);
//}

void benchmark_quicked(
    align_input_t* const align_input, 
    const int window_size, 
    const int overlap_size, 
    const int bandwidth, 
    const int force_scalar, 
    const int hew_threshold, 
    const int hew_percentage) {
  
  quicked_aligner_t aligner;                          // Aligner object
  quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.
  params.external_timer = true;
  params.external_allocator = align_input->mm_allocator;

  params.window_size = window_size;                     
  params.overlap_size = overlap_size;                     
  params.force_scalar = force_scalar;
  params.bandwidth = bandwidth;
  params.hew_threshold[0] = hew_threshold;
  params.hew_threshold[1] = hew_threshold;
  params.hew_percentage[0] = hew_percentage;
  params.hew_percentage[1] = hew_percentage;

  quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  //cigar_t* const cigar = cigar_new(2*max_cigar_length,align_input->mm_allocator);

  aligner.timer = &align_input->timer;
  aligner.timer_windowed_s = &align_input->timer_windowed_s;
  aligner.timer_windowed_l = &align_input->timer_windowed_l;
  aligner.timer_banded = &align_input->timer_banded;
  aligner.timer_align = &align_input->timer_align;
  
  // Align
  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  // DEBUG
  if (align_input->debug_flags) { // TODO
    //benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}

void benchmark_banded(
    align_input_t* const align_input, 
    const int bandwidth, 
    const int only_score) {
  
  quicked_aligner_t aligner;                          // Aligner object
  quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.
  params.external_timer = true;

  params.algo = BANDED;                               // Select the algorithm: Banded
  params.only_score = only_score;                      
  params.bandwidth = bandwidth;                       // Banded needs a bandwidth
  params.external_allocator = align_input->mm_allocator;

  quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters

  aligner.timer = &align_input->timer;
  aligner.timer_windowed_s = &align_input->timer_windowed_s;
  aligner.timer_windowed_l = &align_input->timer_windowed_l;
  aligner.timer_banded = &align_input->timer_banded;
  aligner.timer_align = &align_input->timer_align;
  
  // Align
  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  // DEBUG
  if (align_input->debug_flags) { // TODO
    // benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}

void benchmark_hirschberg(
    align_input_t* const align_input, 
    const int bandwidth) {
  
  quicked_aligner_t aligner;                          // Aligner object
  quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.
  params.external_timer = true;

  params.algo = HIRSCHBERG;                               // Select the algorithm: Banded
  params.bandwidth = bandwidth;                       // Banded needs a bandwidth
  params.external_allocator = align_input->mm_allocator;

  quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters

  aligner.timer = &align_input->timer;
  aligner.timer_windowed_s = &align_input->timer_windowed_s;
  aligner.timer_windowed_l = &align_input->timer_windowed_l;
  aligner.timer_banded = &align_input->timer_banded;
  aligner.timer_align = &align_input->timer_align;
  
  // Align
  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  // DEBUG
  if (align_input->debug_flags) { // TODO
    // benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}


void benchmark_windowed(
    align_input_t* const align_input, 
    const int window_size, 
    const int overlap_size, 
    const int force_scalar, 
    const int only_score) {
  
  quicked_aligner_t aligner;                          // Aligner object
  quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.
  params.external_timer = true;

  params.algo = WINDOWED;                               // Select the algorithm: WindowEd
  params.only_score = only_score;                      
  params.window_size = window_size;                     
  params.overlap_size = overlap_size; 
  params.force_scalar = force_scalar;                     
  params.external_allocator = align_input->mm_allocator;

  quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters

  aligner.timer = &align_input->timer;
  aligner.timer_windowed_s = &align_input->timer_windowed_s;
  aligner.timer_windowed_l = &align_input->timer_windowed_l;
  aligner.timer_banded = &align_input->timer_banded;
  aligner.timer_align = &align_input->timer_align;
  
  // Align
  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  // DEBUG
  if (align_input->debug_flags) { // TODO
    // benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}

