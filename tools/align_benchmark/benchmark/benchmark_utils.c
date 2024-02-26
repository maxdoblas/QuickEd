/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
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
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Benchmark utils
 */

#include "benchmark/benchmark_utils.h"
#include "score_matrix.h"
#include "edit/edit_dp.h"

/*
 * Setup
 */
void benchmark_align_input_clear(
    align_input_t* const align_input) {
  // Alignment form
  align_input->ends_free = false;
  align_input->pattern_begin_free = 0;
  align_input->text_begin_free = 0;
  align_input->pattern_end_free = 0;
  align_input->text_end_free = 0;
  // Output
  align_input->output_file = NULL;
  align_input->output_full = false;
  // Accuracy Stats
  counter_reset(&(align_input->align));
  counter_reset(&(align_input->align_correct));
  counter_reset(&(align_input->align_score));
  counter_reset(&(align_input->align_score_total));
  counter_reset(&(align_input->align_score_diff));
  counter_reset(&(align_input->align_cigar));
  counter_reset(&(align_input->align_bases));
  counter_reset(&(align_input->align_matches));
  counter_reset(&(align_input->align_mismatches));
  counter_reset(&(align_input->align_del));
  counter_reset(&(align_input->align_ins));
}


/*
 * Display
 */
void benchmark_print_alignment(
    FILE* const stream,
    align_input_t* const align_input,
    const int score_computed,
    cigar_t* const cigar_computed,
    const int score_correct,
    cigar_t* const cigar_correct) {
  // Print Sequence
  fprintf(stream,"ALIGNMENT (#%d)\n",align_input->sequence_id);
  fprintf(stream,"  PATTERN  %s\n",align_input->pattern);
  fprintf(stream,"  TEXT     %s\n",align_input->text);
  // Print CIGARS
  if (cigar_computed != NULL && score_computed != -1) {
    fprintf(stream,"    COMPUTED\tscore=%d\t",score_computed);
    cigar_print(stream,cigar_computed,true,align_input->mm_allocator);
    fprintf(stream,"\n");
  }
  if (cigar_computed != NULL) {
    cigar_print_pretty(stream,cigar_computed,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,align_input->mm_allocator);
  }
  if (cigar_correct != NULL && score_correct != -1) {
    fprintf(stream,"    CORRECT \tscore=%d\t",score_correct);
    cigar_print(stream,cigar_correct,true,align_input->mm_allocator);
    fprintf(stream,"\n");
  }
  if (cigar_correct != NULL) {
    cigar_print_pretty(stream,cigar_correct,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,align_input->mm_allocator);
  }
}
void benchmark_print_output_lite(
    FILE* const stream,
    const int score,
    cigar_t* const cigar) {
  // Retrieve CIGAR
  const bool cigar_null = (cigar->begin_offset >= cigar->end_offset);
  char* cigar_str = NULL;
  if (!cigar_null) {
    int cigar_length = 2*(cigar->end_offset-cigar->begin_offset)+10;
    cigar_str = malloc(cigar_length);
    cigar_sprint(cigar_str,cigar_length,cigar,true);
  }
  // Print
  fprintf(stream,"%d\t%s\n",score,(cigar_null) ? "-" : cigar_str);
  // Free
  if (!cigar_null) free(cigar_str);
}
void benchmark_print_output_full(
    FILE* const stream,
    align_input_t* const align_input,
    const int score,
    cigar_t* const cigar) {
  // Retrieve CIGAR
  const bool cigar_null = (cigar->begin_offset >= cigar->end_offset);
  char* cigar_str = NULL;
  if (!cigar_null) {
    int cigar_length = 2*(cigar->end_offset-cigar->begin_offset);
    cigar_str = malloc(2*(cigar->end_offset-cigar->begin_offset));
    cigar_sprint(cigar_str,cigar_length,cigar,true);
  }
  // Print
  fprintf(stream,"%d\t%d\t%d\t%s\t%s\t%s\n",
      align_input->pattern_length,     // Pattern length
      align_input->text_length,        // Text length
      score,                           // Alignment score
      align_input->pattern,            // Pattern sequence
      align_input->text,               // Text sequence
      (cigar_null) ? "-" : cigar_str); // CIGAR
  // Free
  if (!cigar_null) free(cigar_str);
}
void benchmark_print_output(
    align_input_t* const align_input,
    const bool score_only,
    cigar_t* const cigar) {
  if (align_input->output_file) {
    // Compute score
    int score = -1;
    if (score_only) {
      score = cigar->score;
    } else if (cigar->begin_offset < cigar->end_offset) {
      score = cigar_score_edit(cigar);
      cigar->score = score;
    }
    // Print summary
    if (align_input->output_full) {
      benchmark_print_output_full(align_input->output_file,align_input,score,cigar);
    } else {
      benchmark_print_output_lite(align_input->output_file,score,cigar);
    }
  }
}
void quicked_print_output(
    align_input_t* const align_input,
    const bool score_only,
    char* const cigar,
    int score) {
  if (align_input->output_file) {
    // Print summary
    if (align_input->output_full) {
      fprintf(align_input->output_file,"%d\t%d\t%d\t%s\t%s\t%s\n",
              align_input->pattern_length,     // Pattern length
              align_input->text_length,        // Text length
              score,                           // Alignment score
              align_input->pattern,            // Pattern sequence
              align_input->text,               // Text sequence
              (score_only) ? "-" : cigar); // CIGAR
    } else {
      fprintf(align_input->output_file,"%d\t%s\n",score,(score_only) ? "-" : cigar);
    }
  }
}
void quicked_print_output_error(
    align_input_t* const align_input,
    const bool score_only) {
  if (align_input->output_file) {
    // Print summary
    if (align_input->output_full) {
      fprintf(align_input->output_file,"%d\t%d\t%s\t%s\t%s\t%s\n",
              align_input->pattern_length,     // Pattern length
              align_input->text_length,        // Text length
              "ERROR",                         // Alignment score
              align_input->pattern,            // Pattern sequence
              align_input->text,               // Text sequence
              (score_only) ? "-" : "ERROR"); // CIGAR
    } else {
      fprintf(align_input->output_file,"%s\t%s\n","ERROR",(score_only) ? "-" : "ERROR");
    }
  }
}
/*
 * Stats
 */
void benchmark_print_stats(
    FILE* const stream,
    align_input_t* const align_input) {
  // General stats
  fprintf(stream,"[Accuracy]\n");
  fprintf(stream," => Alignments.Correct     ");
  counter_print(stream,&align_input->align_correct,&align_input->align,"alg       ",true);
  // Score check
  if (align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) {
    fprintf(stream," => Score.Correct          ");
    counter_print(stream,&align_input->align_score,&align_input->align,"alg       ",true);
    fprintf(stream,"   => Score.Total          ");
    counter_print(stream,&align_input->align_score_total,NULL,"score uds.",true);
    fprintf(stream,"     => Score.Diff         ");
    counter_print(stream,&align_input->align_score_diff,&align_input->align_score_total,"score uds.",true);
  }
  // Alignment check
  if (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT) {
    fprintf(stream," => CIGAR.Correct          ");
    counter_print(stream,&align_input->align_cigar,&align_input->align,"alg       ",true);
  }
  // CIGAR stats
  fprintf(stream," => CIGAR.Breakdown        \n");
  fprintf(stream,"   => CIGAR.Matches        ");
  counter_print(stream,&align_input->align_matches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Mismatches     ");
  counter_print(stream,&align_input->align_mismatches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Insertions     ");
  counter_print(stream,&align_input->align_ins,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Deletions      ");
  counter_print(stream,&align_input->align_del,&align_input->align_bases,"bases     ",true);
}

