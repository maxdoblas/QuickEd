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

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "alignment/bpm_banded.h"
#include "alignment/bpm.h"
#include "utils/dna_text.h"
#include <immintrin.h>

/*
 * Constants
 */
#define BPM_W64_LENGTH UINT64_LENGTH
#define BPM_W64_SIZE   UINT64_SIZE

void bpm_compute_matrix_hirschberg(
    char* const text,
    char* const text_r,
    const int64_t text_length,
    char* const pattern,
    char* const pattern_r,
    const int64_t pattern_length,
    const uint64_t cutoff_score,
    const cigar_t* cigar_out,
    mm_allocator_t* const mm_allocator){

  const int64_t k_end = ABS(((int64_t)text_length)-(int64_t)(pattern_length))+1;
  const int64_t cutoff_score_real = MAX(MAX(k_end,cutoff_score),65);
  const int64_t sequence_length_diff = pattern_length - text_length;
  const int64_t relative_cutoff_score = DIV_CEIL((cutoff_score_real - ABS(sequence_length_diff)), 2);
  int64_t prolog_column_blocks;
  int64_t effective_bandwidth_blocks;
  if(sequence_length_diff >= 0){
    prolog_column_blocks = DIV_CEIL(relative_cutoff_score, BPM_W64_LENGTH);
    effective_bandwidth_blocks = DIV_CEIL(relative_cutoff_score + sequence_length_diff, BPM_W64_LENGTH) + 1 + prolog_column_blocks;
  } else {
    prolog_column_blocks = DIV_CEIL(relative_cutoff_score - sequence_length_diff, BPM_W64_LENGTH);
    effective_bandwidth_blocks = DIV_CEIL(relative_cutoff_score, BPM_W64_LENGTH) + 1 + prolog_column_blocks;
  }

  int64_t alignment_footprint = effective_bandwidth_blocks*text_length*BPM_W64_SIZE*2;

  if (alignment_footprint > BUFFER_SIZE_16M){ // divide the alignment in 2
    //printf("-------------------------------------------------------------\n");
    //printf("divining sequence of size: %ld\n",alignment_footprint);
    //printf("text, pattern: %ld, %lu\n",text_length, pattern_length);
    //printf("cutoff_score, effective_bandwidth_blocks: %lu, %lu\n",cutoff_score, effective_bandwidth_blocks);
    const int64_t text_len = (text_length+1)/2;
    const int64_t text_len_r = text_length - text_len;

    const int64_t pattern_len = pattern_length;

    banded_pattern_t banded_pattern;
    banded_pattern_t banded_pattern_r;
    banded_pattern_compile(
        &banded_pattern,pattern,
        pattern_length,mm_allocator);
    banded_pattern_compile(
        &banded_pattern_r,pattern_r,
        pattern_length,mm_allocator);

    banded_matrix_t banded_matrix, banded_matrix_r;


    // Compute left side (for getting the central column)
    banded_matrix_allocate_cutoff_score(
        &banded_matrix,pattern_length,
        text_length,cutoff_score,mm_allocator);
    
    banded_compute_cutoff_score(
        &banded_matrix,&banded_pattern,text,
        text_length, text_len, cutoff_score);
    //printf("computing_l\n");
    //printf("text_length, text_len, cutoff_score: %lu, %lu, %lu\n",text_length, text_len, cutoff_score);


    // Compute right side (for getting the central column)
    banded_matrix_allocate_cutoff_score(
        &banded_matrix_r,pattern_length,
        text_length,cutoff_score,mm_allocator);
    //printf("computing_r\n");
    //printf("text_length, text_len_r, cutoff_score: %lu, %lu, %lu\n",text_length, text_len, cutoff_score);


    banded_compute_cutoff_score(
        &banded_matrix_r,&banded_pattern_r,text_r,
        text_length, text_len_r, cutoff_score);

    int64_t first_block_band_pos_v = text_len < prolog_column_blocks*BPM_W64_LENGTH ? 0 : (text_len/BPM_W64_LENGTH) - (prolog_column_blocks);
    int64_t first_block_band_pos_v_r = text_len_r < prolog_column_blocks*BPM_W64_LENGTH ? 0 : (text_len_r/BPM_W64_LENGTH) - (prolog_column_blocks);
    //printf("sequence_length_diff, prolog_column_blocks, first_block_band_pos_v, first_block_band_pos_v_r = %ld, %ld, %ld, %ld\n", sequence_length_diff, prolog_column_blocks, first_block_band_pos_v, first_block_band_pos_v_r);
    //printf(".lower_block, .higher_block, _r.lower_block, _r.higher_block = %ld, %ld, %ld, %ld\n", banded_matrix.lower_block, banded_matrix.higher_block, banded_matrix_r.lower_block, banded_matrix_r.higher_block);

    int64_t bottom_cell;
    int64_t higher_cell, higher_cell_r;
    int64_t starting_pos;
    const int64_t bottom_pos = banded_matrix.lower_block*64 + 63 + first_block_band_pos_v*64;
    const int64_t bottom_pos_r = (pattern_len-1) - (banded_matrix_r.higher_block*64 + 63 + first_block_band_pos_v_r*64);
    const int64_t higher_pos = banded_matrix.higher_block*64 + 63 + first_block_band_pos_v*64;
    const int64_t higher_pos_r = (pattern_len-1) - (banded_matrix_r.lower_block*64 + 63 + first_block_band_pos_v_r*64);

    //printf("bottom_pos, bottom_pos_r, higher_pos, higher_pos_r = %ld, %ld, %ld, %ld\n", bottom_pos, bottom_pos_r, higher_pos, higher_pos_r);
    
    // select lower cell
    if (bottom_pos > bottom_pos_r){
      bottom_cell = banded_matrix.lower_block*64 + 63;
      starting_pos = bottom_pos;
    }else{
      bottom_cell = bottom_pos_r - first_block_band_pos_v*64;
      starting_pos = bottom_pos_r;
    }

    // select higher cell
    if (higher_pos < higher_pos_r){
      higher_cell = banded_matrix.higher_block*64 + 63;
      higher_cell_r = (pattern_len-1) - higher_pos - first_block_band_pos_v_r*64;
    }else{
      higher_cell = higher_pos_r - first_block_band_pos_v*64;
      higher_cell_r = banded_matrix_r.lower_block*64 + 63;
    }
    const int64_t number_of_cells = higher_cell - bottom_cell + 2;
    //printf("bottom_cell, higher_cell, bottom_cell_r, higher_cell_r, starting_pos, number_of_cells : %ld,%ld,%ld,%ld,%ld,%ld \n", bottom_cell, higher_cell, bottom_cell_r, higher_cell_r, starting_pos,number_of_cells);

    int32_t* cell_score = (int32_t*)mm_allocator_malloc(mm_allocator,number_of_cells*sizeof(int32_t));
    int32_t* cell_score_r = (int32_t*)mm_allocator_malloc(mm_allocator,number_of_cells*sizeof(int32_t));

    // compute scores of the left side
    cell_score[0] = 0;
    for (uint64_t i = 0; i < number_of_cells; i++){
      const uint64_t block = (bottom_cell+i)/BPM_W64_LENGTH;
      const uint64_t cell = (bottom_cell+i)%BPM_W64_LENGTH;
      cell_score[i+1] = cell_score[i] + ((banded_matrix.Pv[block] >> cell) & 0x1ULL) - ((banded_matrix.Mv[block] >> cell) & 0x1ULL);
    }
    // compute scores of the right side
    cell_score_r[0] = 0;
    for (uint64_t i = 0; i < number_of_cells; i++){
      const uint64_t block = (higher_cell_r+i)/BPM_W64_LENGTH;
      const uint64_t cell = (higher_cell_r+i)%BPM_W64_LENGTH;
      cell_score_r[i+1] = cell_score_r[i] + ((banded_matrix_r.Pv[block] >> cell) & 0x1ULL) - ((banded_matrix_r.Mv[block] >> cell) & 0x1ULL);
    }

    // search the middle joint cell
    int64_t smaller_pos = 0;
    int64_t smaller_score = cell_score_r[number_of_cells-1] + cell_score[0];
    for (uint64_t i = 1; i < number_of_cells; i++){
      int64_t new_score = cell_score_r[number_of_cells-1-i] + cell_score[i];
      if (new_score < smaller_score){
        smaller_pos = i;
        smaller_score = new_score;
      }
    }
  
    
    int64_t pattern_length_left = starting_pos + smaller_pos;
    int64_t pattern_length_right = pattern_length - pattern_length_left;

    int64_t score_pos_r = DIV_CEIL(pattern_length_right,BPM_W64_LENGTH)*BPM_W64_LENGTH - (higher_cell_r+ first_block_band_pos_v_r*64);
    int64_t score_pos_l = DIV_CEIL(pattern_length_left,BPM_W64_LENGTH)*BPM_W64_LENGTH - (bottom_cell+ first_block_band_pos_v*64);

    char* pattern_r_left = pattern_r + pattern_length_right;
    char* pattern_right = pattern + pattern_length_left;

    int64_t text_length_right = text_length - text_len;
    char* text_right = text + text_len;
    char* text_r_left = text_r + text_length_right;

    int64_t score_r = cell_score_r[number_of_cells-1-smaller_pos] - cell_score_r[score_pos_r] + banded_matrix_r.scores[(pattern_length_right)/BPM_W64_LENGTH];
    int64_t score_l = cell_score[smaller_pos] - cell_score[score_pos_l] + banded_matrix.scores[(pattern_length_left)/BPM_W64_LENGTH];

    // Free
    banded_pattern_free(&banded_pattern,mm_allocator);
    banded_pattern_free(&banded_pattern_r,mm_allocator);
    banded_matrix_free_cutoff(&banded_matrix,mm_allocator);
    banded_matrix_free_cutoff(&banded_matrix_r,mm_allocator);
    mm_allocator_free(mm_allocator,cell_score);
    mm_allocator_free(mm_allocator,cell_score_r);

    // Compute right
    bpm_compute_matrix_hirschberg(
      text_right,
      text_r,
      text_length_right,
      pattern_right,
      pattern_r,
      pattern_length_right,
      score_r,
      cigar_out,
      mm_allocator);

    // Compute left
    bpm_compute_matrix_hirschberg(
      text,
      text_r_left,
      text_len,
      pattern,
      pattern_r_left,
      pattern_length_left,
      score_l,
      cigar_out,
      mm_allocator);

  }else{ // solve the alignment
    //printf("-------------------------------------------------------------\n");
    //printf("Computing sequence of size: %lu\n",alignment_footprint);
    //printf("text, pattern: %lu, %lu\n",text_length, pattern_length);
    //printf("cutoff_score, effective_bandwidth_blocks: %lu, %lu\n",cutoff_score, effective_bandwidth_blocks);

    banded_pattern_t banded_pattern;
    banded_matrix_t banded_matrix;
    
    // Allocate
    banded_pattern_compile(
        &banded_pattern,pattern,
        pattern_length,mm_allocator);
    banded_matrix_allocate_cutoff(
        &banded_matrix,pattern_length,
        text_length,cutoff_score,mm_allocator);
    
    // Align
    banded_compute_cutoff(
        &banded_matrix,&banded_pattern,text,
        text_length,cutoff_score);
    // Merge cigar    
    cigar_insert(banded_matrix.cigar,cigar_out);
    // free variables
    banded_pattern_free(&banded_pattern,mm_allocator);
    banded_matrix_free_cutoff(&banded_matrix,mm_allocator);
  }
}

