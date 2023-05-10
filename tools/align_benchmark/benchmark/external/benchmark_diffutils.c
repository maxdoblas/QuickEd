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

/*
 * Benchmark DiffUtils Myers-O(ND)
 */
#include <limits.h>
#include <stdbool.h>
#include "external/diffutils/minmax.h"

#define ELEMENT char
#define EQUAL(x,y) ((x) == (y))
#define OFFSET int
#define EXTRA_CONTEXT_FIELDS char *deletions; char *insertions;
#define NOTE_DELETE(ctxt,xoff) ctxt->deletions[xoff] = 'D'
#define NOTE_INSERT(ctxt,yoff) ctxt->insertions[yoff] = 'I'
#define USE_HEURISTIC 1
#include "external/diffutils/diffseq.h"

void benchmark_diffutils(
    align_input_t* const align_input,
    const bool find_minimal) {
  // Allocate CIGAR
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t cigar = {
      .operations = malloc(max_cigar_length),
      .begin_offset = 0,
      .end_offset = 0,
  };
  // Prepare diff-alignment
  struct context ctxt;
  ctxt.xvec = align_input->pattern;
  ctxt.yvec = align_input->text;
  const int diags = (align_input->pattern_length + align_input->text_length + 3);
  ctxt.fdiag = malloc(diags * 2 * sizeof(int));
  ctxt.bdiag = ctxt.fdiag + diags;
  ctxt.fdiag += align_input->pattern_length + 1;
  ctxt.bdiag += align_input->pattern_length + 1;
  ctxt.heuristic = !find_minimal;
  ctxt.too_expensive = 1000000; // 4096;
  void* ops_mem = malloc(2*max_cigar_length);
  ctxt.deletions = ops_mem;
  ctxt.insertions = ops_mem + max_cigar_length;
  memset(ops_mem,0,2*max_cigar_length);
  // Align
  timer_start(&align_input->timer);
  compareseq(0,align_input->pattern_length,0,align_input->text_length,find_minimal,&ctxt);
  timer_stop(&align_input->timer);
  // Adapt
  if (align_input->debug_flags || align_input->output_file) {
    // Adapt CIGAR
    int cigar_pos = 0, pattern_pos = 0, text_pos = 0;
    while (pattern_pos < align_input->pattern_length ||
           text_pos < align_input->text_length) {
      if (ctxt.deletions[pattern_pos] == 'D') {
        cigar.operations[cigar_pos++] = 'D';
        pattern_pos++;
      } else if (ctxt.insertions[text_pos] == 'I') {
        cigar.operations[cigar_pos++] = 'I';
        text_pos++;
      } else {
        cigar.operations[cigar_pos++] = 'M';
        pattern_pos++;
        text_pos++;
      }
    }
    cigar.begin_offset = 0;
    cigar.end_offset = cigar_pos;
  }
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,&cigar);
  }
  // Free
  free(ctxt.fdiag - (align_input->pattern_length + 1));
  free(ops_mem);
  free(cigar.operations);
}
