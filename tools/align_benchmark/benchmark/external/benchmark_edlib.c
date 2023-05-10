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
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_check.h"

#include "external/edlib/edlib/include/edlib.h"

/*
 * Benchmark EdLib
 */
void benchmark_edlib_adapt_cigar(
    align_input_t* const align_input,
    char* const edlib_cigar,
    cigar_t* const cigar) {
  // Decode all CIGAR operations
  const int cigar_length = strlen(edlib_cigar);
  int chars_read = 0;
  while (chars_read < cigar_length) {
    // Read operation
    int length, n;
    char operation;
    sscanf(edlib_cigar+chars_read,"%d%c%n",&length,&operation,&n);
    chars_read += n;
    // Adapt operation encoding
    if (operation=='=') operation = 'M';
    else if (operation=='D') operation = 'I';
    else if (operation=='I') operation = 'D';
    // Dump operation
    int i;
    for (i=0;i<length;++i) {
      cigar->operations[(cigar->end_offset)++] = operation;
    }
  }
}
void benchmark_edlib(align_input_t* const align_input) {
  // Parameters
  EdlibAlignResult result;
  char* edlib_cigar = NULL;
  // Align
  timer_start(&align_input->timer);
  result = edlibAlign(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      edlibNewAlignConfig(-1,EDLIB_MODE_NW,EDLIB_TASK_PATH,NULL,0));
  edlib_cigar = edlibAlignmentToCigar(
      result.alignment,result.alignmentLength,EDLIB_CIGAR_EXTENDED); // Traceback
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  cigar.operations = malloc(align_input->pattern_length+align_input->text_length);
  cigar.begin_offset = 0;
  cigar.end_offset = 0;
  benchmark_edlib_adapt_cigar(align_input,edlib_cigar,&cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,&cigar);
  }
  // Free
  free(edlib_cigar);
  free(cigar.operations);
  edlibFreeAlignResult(result);
}
