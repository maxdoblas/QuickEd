/*
 *                             The MIT License
 *
 * This file is part of QuickEd library.
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
#include "edit/edit_bpm.h"
#include "../../../quicked/quicked.h"
#include "quicked_utils/include/commons.h"

#include "external/edlib/edlib/include/edlib.h"
#include "external/astar-pairwise-aligner/astarpa.h"
#include "external/SneakySnake/SneakySnake/SneakySnake.h"
/*
 * Benchmark Edit
 */
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

  aligner.timer = &align_input->timer;
  aligner.timer_windowed_s = &align_input->timer_windowed_s;
  aligner.timer_windowed_l = &align_input->timer_windowed_l;
  aligner.timer_banded = &align_input->timer_banded;
  aligner.timer_align = &align_input->timer_align;
  
  // Align
  quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);
  
  // DEBUG
  if (align_input->debug_flags) {
    cigar2_t* const cigar = cigar_new_2(
      align_input->pattern_length + align_input->text_length,
      align_input->mm_allocator
    );
    cigar_to_operations_2(cigar, aligner.cigar, strlen(aligner.cigar));

    benchmark_check_alignment(align_input,cigar);

    cigar_free_2(cigar, align_input->mm_allocator);
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
  if (align_input->debug_flags) {
    cigar2_t* const cigar = cigar_new_2(
      align_input->pattern_length + align_input->text_length,
      align_input->mm_allocator
    );
    cigar_to_operations_2(cigar, aligner.cigar, strlen(aligner.cigar));

    benchmark_check_alignment(align_input,cigar);

    cigar_free_2(cigar, align_input->mm_allocator);
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
  quicked_status_t status = quicked_align(&aligner, align_input->pattern, align_input->pattern_length, align_input->text, align_input->text_length);

   if (quicked_check_error(status)){
    fprintf(stderr, "%s", quicked_status_msg(status));
   }
  
  // DEBUG
  if (align_input->debug_flags) {
    cigar2_t* const cigar = cigar_new_2(
      align_input->pattern_length + align_input->text_length,
      align_input->mm_allocator
    );
    cigar_to_operations_2(cigar, aligner.cigar, strlen(aligner.cigar));

    benchmark_check_alignment(align_input,cigar);

    cigar_free_2(cigar, align_input->mm_allocator);
  }
  // Output
  if (align_input->output_file) {
    if (quicked_check_error(status)){
      quicked_print_output_error(align_input,false);
    }else{
      quicked_print_output(align_input,false,aligner.cigar,aligner.score);
    }
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
  if (align_input->debug_flags) {
    cigar2_t* const cigar = cigar_new_2(
      align_input->pattern_length + align_input->text_length,
      align_input->mm_allocator
    );
    cigar_to_operations_2(cigar, aligner.cigar, strlen(aligner.cigar));

    benchmark_check_alignment(align_input,cigar);

    cigar_free_2(cigar, align_input->mm_allocator);
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,aligner.cigar,aligner.score);
  }
  // Free
  quicked_free(&aligner);                 // Free whatever memory the aligner allocated

}


void benchmark_edit_bpm(
    align_input_t* const align_input) {
  // Allocate
  bpm_pattern_t bpm_pattern;
  edit_bpm_pattern_compile(
      &bpm_pattern,align_input->pattern,
      align_input->pattern_length,align_input->mm_allocator);
  bpm_matrix_t bpm_matrix;
  edit_bpm_matrix_allocate(
      &bpm_matrix,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  edit_bpm_compute(
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
  edit_bpm_pattern_free(&bpm_pattern,align_input->mm_allocator);
  edit_bpm_matrix_free(&bpm_matrix,align_input->mm_allocator);
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
  cigar2_t* const cigar = cigar_new_2(
      pattern_length+text_length, align_input->mm_allocator);
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
  cigar_free_2(cigar,align_input->mm_allocator);
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
  cigar2_t* const cigar = cigar_new_2(
      pattern_length+text_length,align_input->mm_allocator);
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
  cigar_free_2(cigar,align_input->mm_allocator);
}


// External

void benchmark_edlib(align_input_t* const align_input,
    const int bandwidth) {
  
  const int pattern_length = align_input->pattern_length;
  const int text_length = align_input->text_length;
  int64_t bandwidth_k = bandwidth != -1 ? (MAX(text_length,pattern_length)*bandwidth)/100 : -1;
  // Parameters
  EdlibAlignResult result;
  char* edlib_cigar = NULL;
  // Align
  timer_start(&align_input->timer);
  result = edlibAlign(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      edlibNewAlignConfig(bandwidth_k,EDLIB_MODE_NW,EDLIB_TASK_PATH,NULL,0));
  edlib_cigar = edlibAlignmentToCigar(
      result.alignment,result.alignmentLength,EDLIB_CIGAR_EXTENDED); // Traceback
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  uint64_t cigar_len = strlen(edlib_cigar);
  for(uint64_t i = 0; i < cigar_len;i++){
    char operation = edlib_cigar[i];
    if (operation=='=') edlib_cigar[i] = 'M';
    else if (operation=='D') edlib_cigar[i] = 'I';
    else if (operation=='I') edlib_cigar[i] = 'D';
  }
  // Output
  if (align_input->output_file) {
    quicked_print_output(align_input,false,edlib_cigar,result.editDistance);
  }
  // Free
  free(edlib_cigar);
  edlibFreeAlignResult(result);
}


/*
void benchmark_scrooge(align_input_t* const align_input) {
  // Parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t* const cigar = cigar_new(2*max_cigar_length);
  // Align
  uint64_t time_ns = 0;
  benchmark_scrooge_bridge(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      cigar->operations,&cigar->end_offset,&time_ns);
  counter_add(&align_input->timer.time_ns,time_ns);
  cigar_score_edit(cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,false,cigar);
  }
  // Free
  cigar_free(cigar);
}*/

void benchmark_astarpa(align_input_t* const align_input)
{
  const int pattern_length = align_input->pattern_length;
  const int text_length    = align_input->text_length;
  uint64_t cigar_len; 
  uint8_t* cigar; 
  timer_start(&align_input->timer);
	uint64_t score = astarpa((const uint8_t*)align_input->pattern, pattern_length, 
                      (const uint8_t*)align_input->text,    text_length,    &cigar, &cigar_len);
  timer_stop(&align_input->timer);
  //adapt cigar
	for(uint64_t i = 0; i < cigar_len;i++){
    char operation = cigar[i];
    if (operation=='=') cigar[i] = 'M';
    else if (operation=='D') cigar[i] = 'I';
    else if (operation=='I') cigar[i] = 'D';
  }/*
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }*/
  if (align_input->output_file) {
    quicked_print_output(align_input,false,cigar,score);
  }
  astarpa_free_cigar(cigar);
}

void benchmark_astarpa2_simple(align_input_t* const align_input)
{
  const int pattern_length = align_input->pattern_length;
  const int text_length    = align_input->text_length;
  uint64_t cigar_len; 
  uint8_t* cigar; 
  timer_start(&align_input->timer);
	uint64_t score = astarpa2_simple((const uint8_t*)align_input->pattern, pattern_length, 
                                   (const uint8_t*)align_input->text,    text_length,    &cigar, &cigar_len);
  timer_stop(&align_input->timer);
  //adapt cigar
	for(uint64_t i = 0; i < cigar_len;i++){
    char operation = cigar[i];
    if (operation=='=') cigar[i] = 'M';
    else if (operation=='D') cigar[i] = 'I';
    else if (operation=='I') cigar[i] = 'D';
  }
  /*
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }*/
  if (align_input->output_file) {
    quicked_print_output(align_input,false,cigar,score);
  }
  astarpa_free_cigar(cigar);
}

void benchmark_astarpa2_full(align_input_t* const align_input)
{
  const int pattern_length = align_input->pattern_length;
  const int text_length    = align_input->text_length;
  uint64_t cigar_len; 
  uint8_t* cigar; 
  timer_start(&align_input->timer);
	uint64_t score = astarpa2_full((const uint8_t*)align_input->pattern, pattern_length, 
                                 (const uint8_t*)align_input->text,    text_length,    &cigar, &cigar_len);
  timer_stop(&align_input->timer);
  //adapt cigar
	for(uint64_t i = 0; i < cigar_len;i++){
    char operation = cigar[i];
    if (operation=='=') cigar[i] = 'M';
    else if (operation=='D') cigar[i] = 'I';
    else if (operation=='I') cigar[i] = 'D';
  }/*
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }*/
  if (align_input->output_file) {
    quicked_print_output(align_input,false,cigar,score);
  }
  astarpa_free_cigar(cigar);
}


void benchmark_sneakysnake(align_input_t* const align_input)
{
  const int pattern_length  = align_input->pattern_length;
  const int text_length     = align_input->text_length;
  const int min_length      = (pattern_length < text_length) ? pattern_length : text_length; 
  const int iteration_num   = (0.21f * (float)min_length); 
  timer_start(&align_input->timer);
  int score = SneakySnake(min_length, align_input->pattern, align_input->text, min_length, iteration_num, 0, iteration_num);
  timer_stop(&align_input->timer);
  if (align_input->output_file) {
    quicked_print_output(align_input,false, "",score);
  }
}


void benchmark_wavefront(align_input_t* const align_input)
{
  wavefront_aligner_t* const wf_aligner = align_input->wf_aligner;
  timer_start(&align_input->timer);
  wavefront_align(wf_aligner,
                 align_input->pattern, align_input->pattern_length, 
                 align_input->text, align_input->text_length); 
  timer_stop(&align_input->timer);
  int score = align_input->wf_aligner->cigar->score;
  //cigar_t to cigar2_t
  cigar2_t* cigar = (cigar2_t*)malloc(sizeof(cigar2_t));
  cigar->score          = wf_aligner->cigar->score; 
  cigar->end_h          = wf_aligner->cigar->end_h;
  cigar->end_offset     = wf_aligner->cigar->end_offset;
  cigar->end_v          = wf_aligner->cigar->end_v;
  cigar->begin_offset   = wf_aligner->cigar->begin_offset; 
  cigar->max_operations = wf_aligner->cigar->max_operations; 
  cigar->cigar_length   = wf_aligner->cigar->cigar_length; 
  cigar->cigar_buffer   = wf_aligner->cigar->cigar_buffer;
  cigar->operations     = wf_aligner->cigar->operations; 

  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input, cigar);
  }

  if (align_input->output_file) {
    quicked_print_output(align_input, false, "", score);
  }
  free(cigar);
}



