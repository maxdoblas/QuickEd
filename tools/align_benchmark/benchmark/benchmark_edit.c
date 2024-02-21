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
 * Adapt CIGAR
 */
void adapt_cigar_quicked(
    const int pattern_length,
    const int text_length,
    char* const scooge_cigar,
    cigar_t* const cigar) {
  // Init
  int num_edit_operations = 0;
  // Decode all CIGAR operations
  const int cigar_length = strlen(scooge_cigar);
  int chars_read = 0, plen = 0, tlen = 0;
  int length, n, i;
  char operation;
  while (chars_read < cigar_length) {
    // Read operation
    sscanf(scooge_cigar+chars_read,"%d%c%n",&length,&operation,&n);
    chars_read += n;
    // Adapt operation encoding
    if (operation=='=' || operation=='M') {
      operation = 'M'; plen+=length; tlen+=length;
    } else if (operation=='X') {
      operation = 'X'; plen+=length; tlen+=length;
    } else if (operation=='D') {
      operation = 'I'; tlen+=length;
    } else if (operation=='I') {
      operation = 'D'; plen+=length;
    } else {
      fprintf(stderr,"[Scrooge] Computing CIGAR score: Unknown operation '%c'\n",operation);
      exit(1);
    }
    // Dump operation
    for (i=0;i<length;++i) {
      cigar->operations[num_edit_operations++] = operation;
    }
  }
  // Add final indel
  if (plen < pattern_length) {
    for (i=0;i<pattern_length-plen;++i) {
      cigar->operations[num_edit_operations++] = 'D';
    }
  }
  if (tlen < text_length) {
    for (i=0;i<text_length-tlen;++i) {
      cigar->operations[num_edit_operations++] = 'I';
    }
  }
}

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
    align_input_t* const align_input) {
  
  quicked_aligner_t aligner;                          // Aligner object
  quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.

  quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t* const cigar = cigar_new(2*max_cigar_length,align_input->mm_allocator);
  
  // Align
  timer_start(&align_input->timer);

  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) { // TODO
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}

