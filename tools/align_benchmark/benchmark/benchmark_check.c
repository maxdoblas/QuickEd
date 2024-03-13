/*
 *                             The MIT License
 *
 * This file is part of QuickEd library.
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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

#include "benchmark_check.h"
#include "score_matrix.h"
#include "external/edlib/edlib/include/edlib.h"

/*
 * Checker (given the correct CIGAR+SCORE)
 */
void benchmark_check_alignment_using_solution(
    align_input_t* const align_input,
    cigar_t* const cigar_computed,
    const int score_computed,
    cigar_t* const cigar_correct,
    const int score_correct) {
  counter_add(&(align_input->align),1);
  counter_add(&(align_input->align_score_total),ABS(score_computed));
  // Debug
  if (align_input->debug_flags) {
    // Display info
    if (align_input->debug_flags & ALIGN_DEBUG_DISPLAY_INFO) {
      benchmark_print_alignment(stderr,align_input,score_computed,cigar_computed,-1,NULL);
    }
    // Check correct
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_CORRECT) {
      bool correct = cigar_check_alignment(stderr,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length,
          cigar_computed,align_input->verbose);
      if (!correct) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INCORRECT ALIGNMENT\n");
          benchmark_print_alignment(stderr,align_input,-1,cigar_computed,-1,NULL);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_correct),1);
      }
      // CIGAR Stats
      int i;
      counter_add(&(align_input->align_bases),align_input->pattern_length);
      for (i=cigar_computed->begin_offset;i<cigar_computed->end_offset;++i) {
        switch (cigar_computed->operations[i]) {
          case 'M': counter_add(&(align_input->align_matches),1); break;
          case 'X': counter_add(&(align_input->align_mismatches),1); break;
          case 'I': counter_add(&(align_input->align_ins),1); break;
          case 'D': default: counter_add(&(align_input->align_del),1); break;
        }
      }
    }
    // Check score
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) {
      if (score_computed != score_correct) {
        // Print
        if (align_input->verbose) {
          benchmark_print_alignment(
              stderr,align_input,
              score_computed,cigar_computed,
              score_correct,cigar_correct);
          fprintf(stderr,"(#%d)\t INACCURATE SCORE computed=%d\tcorrect=%d\n",
              align_input->sequence_id,score_computed,score_correct);
        }
        counter_add(&(align_input->align_score_diff),ABS(score_computed-score_correct));
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_score),1);
      }
    }
    // Check alignment
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT) {
      if (cigar_cmp(cigar_computed,cigar_correct) != 0) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INACCURATE ALIGNMENT\n");
          benchmark_print_alignment(
              stderr,align_input,
              -1,cigar_computed,
              -1,cigar_correct);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_cigar),1);
      }
    }
  }
}
/*
 * Specialized Checkers
 */
void benchmark_check_alignment_edit(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  EdlibAlignResult result;
  char* edlib_cigar = NULL;

  cigar_t* const cigar = cigar_new(
    align_input->pattern_length + align_input->text_length,
    align_input->mm_allocator
  );

  result = edlibAlign(
    align_input->pattern, align_input->pattern_length,
    align_input->text, align_input->text_length,
    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0)
  );
  edlib_cigar = edlibAlignmentToCigar(
    result.alignment,
    result.alignmentLength,
    EDLIB_CIGAR_EXTENDED
  );

  uint64_t cigar_len = strlen(edlib_cigar);
  for(uint64_t i = 0; i < cigar_len;i++){
    char operation = edlib_cigar[i];
    if (operation=='=') edlib_cigar[i] = 'M';
    else if (operation=='D') edlib_cigar[i] = 'I';
    else if (operation=='I') edlib_cigar[i] = 'D';
  }
  cigar_to_operations(cigar, edlib_cigar, cigar_len);

  const int score_correct = cigar_score_edit(cigar);
  const int score_computed = cigar_score_edit(cigar_computed);

  benchmark_check_alignment_using_solution(
    align_input, cigar_computed, score_computed,
    cigar, score_correct);

  cigar_free(cigar, align_input->mm_allocator);
  free(edlib_cigar);
  edlibFreeAlignResult(result);
}
/*
 * Check
 */
void benchmark_check_alignment(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  // Compute correct CIGAR
  if ((align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) ||
      (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT)) {
    // Check edit
    benchmark_check_alignment_edit(align_input,cigar_computed);
  } else if ((align_input->debug_flags & ALIGN_DEBUG_CHECK_CORRECT) ||
             (align_input->debug_flags & ALIGN_DEBUG_DISPLAY_INFO)) { // Check CORRECT or just DISPLAY
    // Delegate check alignment
    benchmark_check_alignment_using_solution(
        align_input,cigar_computed,-1,NULL,-1);
  }
}

