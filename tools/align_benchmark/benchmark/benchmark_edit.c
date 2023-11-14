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
#include "../../../alignment/bpm.h"
#include "../../../alignment/bpm_banded.h"
#include "../../../alignment/bpm_windowed.h"
#include "../../../alignment/bpm_hirschberg.h"
#include "utils/commons.h"

#define BPM_W64_LENGTH UINT64_LENGTH


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
  
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  const int bandwidth_k = (MAX(text_length,pattern_length)*bandwidth)/100;

  // Allocate
  banded_pattern_t banded_pattern;
  banded_pattern_compile(
      &banded_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  banded_matrix_t banded_matrix;
  banded_matrix_allocate(
      &banded_matrix,align_input->pattern_length,
      align_input->text_length,bandwidth_k,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  banded_compute(
      &banded_matrix,&banded_pattern,align_input->text,
      align_input->text_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,banded_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,banded_matrix.cigar);
  }
  // Free
  banded_pattern_free(&banded_pattern,align_input->mm_allocator);
  banded_matrix_free(&banded_matrix,align_input->mm_allocator);
}

void benchmark_edit_bpm_banded_unaligned(
    align_input_t* const align_input,
    const int bandwidth) {

  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  const int bandwidth_k = (MAX(text_length,pattern_length)*bandwidth)/100;

  // Allocate
  banded_pattern_t banded_pattern;
  banded_pattern_compile(
      &banded_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  banded_matrix_t banded_matrix;
  banded_matrix_allocate_unaligned(
      &banded_matrix,align_input->pattern_length,
      align_input->text_length,bandwidth_k,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  banded_compute_unaligned(
      &banded_matrix,&banded_pattern,align_input->text,
      align_input->text_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,banded_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,banded_matrix.cigar);
  }
  // Free
  banded_pattern_free(&banded_pattern,align_input->mm_allocator);
  banded_matrix_free_unaligned(&banded_matrix,align_input->mm_allocator);
}

void benchmark_edit_bpm_banded_blocking(
    align_input_t* const align_input,
    const int bandwidth) {
  
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  const int bandwidth_k = (MAX(text_length,pattern_length)*bandwidth)/100;

  // Allocate
  banded_pattern_t banded_pattern;
  banded_pattern_compile(
      &banded_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  banded_matrix_t banded_matrix;
  banded_matrix_allocate_blocking(
      &banded_matrix,align_input->pattern_length,
      align_input->text_length,bandwidth_k,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  banded_compute_blocking(
      &banded_matrix,&banded_pattern,align_input->text,
      align_input->text_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,banded_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,banded_matrix.cigar);
  }
  // Free
  banded_pattern_free(&banded_pattern,align_input->mm_allocator);
  banded_matrix_free_blocking(&banded_matrix,align_input->mm_allocator);
}

void benchmark_edit_bpm_banded_cutoff(
    align_input_t* const align_input,
    const int bandwidth) {
  
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  const int bandwidth_k = (MAX(text_length,pattern_length)*bandwidth)/100;

  // Allocate
  banded_pattern_t banded_pattern;
  banded_pattern_compile(
      &banded_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  banded_matrix_t banded_matrix;
  banded_matrix_allocate_cutoff(
      &banded_matrix,align_input->pattern_length,
      align_input->text_length,bandwidth_k,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  banded_compute_cutoff(
      &banded_matrix,&banded_pattern,align_input->text,
      align_input->text_length,bandwidth_k);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,banded_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,banded_matrix.cigar);
  }
  // Free
  banded_pattern_free(&banded_pattern,align_input->mm_allocator);
  banded_matrix_free_cutoff(&banded_matrix,align_input->mm_allocator);
}

void benchmark_edit_bpm_banded_cutoff_score(
    align_input_t* const align_input,
    const int bandwidth) {
  
  const int64_t pattern_length = align_input->pattern_length;
  const int64_t text_length = align_input->text_length;
  const int64_t score = (MAX(text_length,pattern_length)*bandwidth)/100;

  // Allocate
  banded_pattern_t banded_pattern;
  banded_pattern_compile(
      &banded_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  banded_matrix_t banded_matrix;
  banded_matrix_allocate_cutoff_score(
      &banded_matrix,align_input->pattern_length,
      align_input->text_length,score,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  banded_compute_cutoff_score(
      &banded_matrix,&banded_pattern,align_input->text,
      align_input->text_length, align_input->text_length, score);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,banded_matrix.cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,true,banded_matrix.cigar);
  }
  // Free
  banded_pattern_free(&banded_pattern,align_input->mm_allocator);
  banded_matrix_free_cutoff(&banded_matrix,align_input->mm_allocator);
}

void benchmark_edit_bpm_quicked(
    align_input_t* const align_input) {
  
  const int64_t text_len = align_input->text_length;
  const int64_t pattern_len = align_input->pattern_length;

  char* text_r = (char*)mm_allocator_malloc(align_input->mm_allocator,align_input->text_length);
  char* pattern_r = (char*)mm_allocator_malloc(align_input->mm_allocator,align_input->pattern_length);

  reverse_string(align_input->text,text_r,align_input->text_length);
  reverse_string(align_input->pattern,pattern_r,align_input->pattern_length);

  windowed_pattern_t windowed_pattern;
  windowed_matrix_t windowed_matrix;

  timer_start(&align_input->timer);  
  timer_start(&align_input->timer_window_sse);

  windowed_pattern_compile(
      &windowed_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  windowed_matrix_allocate(
      &windowed_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator,
      2);
  // Align
  windowed_compute(
      &windowed_matrix,&windowed_pattern,align_input->text,
      align_input->text_length,align_input->pattern_length,
      2, 1, WINDOW_QUICKED);

  timer_stop(&align_input->timer_window_sse);
  
  //printf("----------------------------------------------------\n");
  //printf("text_len, pattern_len, lim_score = %ld, %ld, %ld\n", text_len,pattern_len,MAX(text_len,pattern_len)/4);
  
  //printf("score wind 2x1 = %ld\n", windowed_matrix.cigar->score);
  int64_t score = windowed_matrix.cigar->score;

  windowed_pattern_free(&windowed_pattern,align_input->mm_allocator);
  windowed_matrix_free(&windowed_matrix,align_input->mm_allocator);

  if(score > MAX(text_len,pattern_len)/4){

    timer_start(&align_input->timer_window_6x2);
    
    windowed_pattern_compile(
        &windowed_pattern,align_input->pattern,
        align_input->pattern_length,align_input->mm_allocator);
    windowed_matrix_allocate(
        &windowed_matrix,align_input->pattern_length,
        align_input->text_length,align_input->mm_allocator,
        6);
    // Align
    windowed_compute(
        &windowed_matrix,&windowed_pattern,align_input->text,
        align_input->text_length,align_input->pattern_length,
        6, 2, WINDOW_SCORE);
    
    align_input->seq_with_6x2 = true;

    //printf("score wind 6x2 = %ld\n", windowed_matrix.cigar->score);
    score = windowed_matrix.cigar->score;
    windowed_pattern_free(&windowed_pattern,align_input->mm_allocator);
    windowed_matrix_free(&windowed_matrix,align_input->mm_allocator);

    windowed_pattern_compile(
        &windowed_pattern,pattern_r,
        align_input->pattern_length,align_input->mm_allocator);
    windowed_matrix_allocate(
        &windowed_matrix,align_input->pattern_length,
        align_input->text_length,align_input->mm_allocator,
        6);
    // Align
    windowed_compute(
        &windowed_matrix,&windowed_pattern,text_r,
        align_input->text_length,align_input->pattern_length,
        6, 2, WINDOW_SCORE);

    align_input->seq_with_6x2_r = true;

    //printf("score wind 6x2 reverse = %ld\n", windowed_matrix.cigar->score);
    score = MIN(score,windowed_matrix.cigar->score);

    windowed_pattern_free(&windowed_pattern,align_input->mm_allocator);
    windowed_matrix_free(&windowed_matrix,align_input->mm_allocator);

    timer_stop(&align_input->timer_window_6x2);

    //printf("score, 0.25 = %ld,%ld\n",score,MAX(text_len,pattern_len)/4);
    if(score > MAX(text_len,pattern_len)/4){

      timer_start(&align_input->timer_banded_15);
      
      banded_pattern_t banded_pattern;
      banded_matrix_t banded_matrix_score;
      // Allocate
      banded_pattern_compile(
        &banded_pattern,align_input->pattern,
        align_input->pattern_length,align_input->mm_allocator);

      score = MAX(text_len,pattern_len)*3/20;

      banded_matrix_allocate_cutoff_score(
        &banded_matrix_score,align_input->pattern_length,
        align_input->text_length,score,
        align_input->mm_allocator);

      banded_compute_cutoff_score(
          &banded_matrix_score,&banded_pattern,align_input->text,
          align_input->text_length, align_input->text_length, score);

      align_input->seqs_with_15 = true;

      int64_t new_score = banded_matrix_score.cigar->score;
      //printf("new score, cutoff = %ld,%ld\n",new_score,score);
      //printf("score band 15 = %ld\n", new_score);


      banded_matrix_free_cutoff(&banded_matrix_score,align_input->mm_allocator);

      timer_stop(&align_input->timer_banded_15);

      while((new_score > MAX(text_len,pattern_len)/4 && score < new_score) || new_score<0){
        score *= 2; 
        timer_start(&align_input->timer_banded_30);

        banded_matrix_allocate_cutoff_score(
          &banded_matrix_score,align_input->pattern_length,
          align_input->text_length,score,
          align_input->mm_allocator); 

        banded_compute_cutoff_score(
            &banded_matrix_score,&banded_pattern,align_input->text,
            align_input->text_length, align_input->text_length, score);
        
        align_input->seqs_with_30 = true;

        //printf("new score, cutoff = %ld,%ld\n",new_score,score);
        //printf("score band 30 = %ld\n", new_score);

        new_score = banded_matrix_score.cigar->score;
        banded_matrix_free_cutoff(&banded_matrix_score,align_input->mm_allocator);

        timer_stop(&align_input->timer_banded_30);

      }

      score = new_score; 
    }

  }

  //printf("Alignment score = %ld\n",score);
  timer_start(&align_input->timer_banded_hirschberg);

  cigar_t cigar_out;
  cigar_out.operations = malloc(align_input->pattern_length+align_input->text_length);
  cigar_out.begin_offset = align_input->pattern_length+align_input->text_length;
  cigar_out.end_offset = align_input->pattern_length+align_input->text_length;

  bpm_compute_matrix_hirschberg(
    align_input->text,
    text_r,
    align_input->text_length,
    align_input->pattern,
    pattern_r,
    align_input->pattern_length,
    score,
    &cigar_out,
    align_input->mm_allocator);

  timer_stop(&align_input->timer_banded_hirschberg);

  timer_stop(&align_input->timer);

  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar_out);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,&cigar_out);
    align_input->diff_scores = (float)(score-cigar_out.score)/(float)(MAX(align_input->text_length,align_input->pattern_length));
  }
  // Free
  free(cigar_out.operations);
}

void benchmark_edit_bpm_windowed(
    align_input_t* const align_input,
    const int window_size,
    const int overlap_size,
    const window_config_t window_config) {
  // Allocate
  windowed_pattern_t windowed_pattern;
  windowed_pattern_compile(
      &windowed_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  windowed_matrix_t windowed_matrix;
  windowed_matrix_allocate(
      &windowed_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator,
      window_size);
  // Align
  timer_start(&align_input->timer);
  windowed_compute(
      &windowed_matrix,&windowed_pattern,align_input->text,
      align_input->text_length,align_input->pattern_length,
      window_size, overlap_size, window_config);
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
  const int bandwidth_k = (MAX(text_length,pattern_length)*bandwidth)/100;
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
      bandwidth_k,cigar);
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

void benchmark_edit_bpm_band_hirschberg(
    align_input_t* const align_input,
    const int bandwidth) {
  // Allocate

  char* text_r = (char*)mm_allocator_malloc(align_input->mm_allocator,align_input->text_length);
  char* pattern_r = (char*)mm_allocator_malloc(align_input->mm_allocator,align_input->pattern_length);

  reverse_string(align_input->text,text_r,align_input->text_length);
  reverse_string(align_input->pattern,pattern_r,align_input->pattern_length);

  // Align
  int score = (MIN(align_input->text_length,align_input->pattern_length)*bandwidth)/100;

  cigar_t cigar_out;
  cigar_out.operations = malloc(align_input->pattern_length+align_input->text_length);
  cigar_out.begin_offset = align_input->pattern_length+align_input->text_length;
  cigar_out.end_offset = align_input->pattern_length+align_input->text_length;

  timer_continue(&align_input->timer);

  bpm_compute_matrix_hirschberg(
    align_input->text,
    text_r,
    align_input->text_length,
    align_input->pattern,
    pattern_r,
    align_input->pattern_length,
    score,
    &cigar_out,
    align_input->mm_allocator);

  timer_stop(&align_input->timer);

  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar_out);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,&cigar_out);
  }
  // Free
  free(cigar_out.operations);
}