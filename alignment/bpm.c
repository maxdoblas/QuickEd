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

#include "alignment/bpm.h"

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "utils/dna_text.h"
#include "common.h"

/*
 * Setup
 */
void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    char* const pattern,
    const int pattern_length,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  const uint64_t PEQ_length = pattern_num_words64*BPM_W64_LENGTH;
  const uint64_t pattern_mod = pattern_length%BPM_W64_LENGTH;
  // Init fields
  bpm_pattern->pattern = pattern;
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words64;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words64*BPM_W64_SIZE;
  const uint64_t PEQ_size = BPM_ALPHABET_LENGTH*aux_vector_size;
  const uint64_t score_size = pattern_num_words64*UINT64_SIZE;
  const uint64_t total_memory = PEQ_size + 3*aux_vector_size + 2*score_size + (pattern_num_words64+1)*UINT64_SIZE;
  void* memory = mm_allocator_malloc(mm_allocator,total_memory);
  bpm_pattern->PEQ = memory; memory += PEQ_size;
  bpm_pattern->P = memory; memory += aux_vector_size;
  bpm_pattern->M = memory; memory += aux_vector_size;
  bpm_pattern->level_mask = memory; memory += aux_vector_size;
  bpm_pattern->score = memory; memory += score_size;
  bpm_pattern->init_score = memory; memory += score_size;
  bpm_pattern->pattern_left = memory;
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = dna_encode(pattern[i]);
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    uint64_t j;
    for (j=0;j<BPM_ALPHABET_LENGTH;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words64-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BPM_W64_MASK;
    bpm_pattern->init_score[i] = BPM_W64_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  for (;i<=pattern_num_words64;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  if (pattern_mod > 0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BPM_W64_MASK;
    bpm_pattern->init_score[top] = BPM_W64_LENGTH;
  }
}
void bpm_pattern_free(
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,bpm_pattern->PEQ);
}
void bpm_matrix_allocate(
    bpm_matrix_t* const bpm_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  bpm_matrix->Mv = Mv;
  bpm_matrix->Pv = Pv;
  // CIGAR
  bpm_matrix->cigar = cigar_new(pattern_length+text_length);
  bpm_matrix->cigar->end_offset = pattern_length+text_length;
}
void bpm_matrix_free(
    bpm_matrix_t* const bpm_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,bpm_matrix->Mv);
  mm_allocator_free(mm_allocator,bpm_matrix->Pv);
  // CIGAR
  cigar_free(bpm_matrix->cigar);
}
/*
 * Edit distance computation using BPM
 */
void bpm_reset_search_cutoff(
    uint8_t* const top_level,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score,
    const uint64_t max_distance) {
  // Calculate the top level (maximum bit-word for cut-off purposes)
  const uint8_t y = (max_distance>0) ? (max_distance+(BPM_W64_LENGTH-1))/BPM_W64_LENGTH : 1;
  *top_level = y;
  // Reset score,P,M
  uint64_t i;
  P[0]=BPM_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<y;++i) {
    P[i]=BPM_W64_ONES;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }
}
void bpm_compute_matrix(
    bpm_matrix_t* const bpm_matrix,
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const uint64_t text_length,
    uint64_t max_distance) {
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  uint64_t* const Pv = bpm_matrix->Pv;
  uint64_t* const Mv = bpm_matrix->Mv;
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level;
  bpm_reset_search_cutoff(&top_level,Pv,Mv,score,init_score,max_distance);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    for (i=0;i<top_level;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];
      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }
    // Cut-off
    const uint8_t last = top_level-1;
    if (score[last]<=max_distance__1 && last<top) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)];
      if (last_score<=max_distance && (MHin || (Peq & 1))) {
        // Init block V
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,top_level);
        const uint64_t next_bdp_idx = bdp_idx+num_words64;
        uint64_t Pv_in = BPM_W64_ONES;
        uint64_t Mv_in = 0;
        Pv[bdp_idx] = BPM_W64_ONES;
        Mv[bdp_idx] = 0;
        const uint64_t mask = level_mask[top_level];
        /* Compute Block */
        BPM_ADVANCE_BLOCK(Peq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
        /* Save Block Pv,Mv */
        Pv[next_bdp_idx]=Pv_in;
        Mv[next_bdp_idx]=Mv_in;
        /* Set score & increment the top level block */
        score[top_level] = last_score + init_score[top_level] + (PHout-MHout);
        ++top_level;
      } else {
        while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
          --top_level;
        }
      }
    } else {
      while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
        --top_level;
      }
    }
  }
  // Return optimal column/distance
  // Check match
  const int64_t current_score = score[top_level-1];
  if (top_level==num_words64 && current_score<=max_distance) {
    bpm_matrix->min_score = score[top_level-1];
    bpm_matrix->min_score_column = text_length-1;
  } else {
    bpm_matrix->min_score = UINT64_MAX;
    bpm_matrix->min_score_column = UINT64_MAX;
  }
}
void bpm_backtrace_matrix(
    bpm_matrix_t* const bpm_matrix,
    const bpm_pattern_t* const bpm_pattern,
    char* const text) {
  // Parameters
  char* const pattern = bpm_pattern->pattern;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t* const Pv = bpm_matrix->Pv;
  const uint64_t* const Mv = bpm_matrix->Mv;
  char* const operations = bpm_matrix->cigar->operations;
  int op_sentinel = bpm_matrix->cigar->end_offset-1;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  int64_t h = bpm_matrix->min_score_column;
  int64_t v = pattern_length - 1;
  while (v >= 0 && h >= 0) {
    const uint8_t block = v / UINT64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h+1,num_words64,block);
    const uint64_t mask = 1UL << (v % UINT64_LENGTH);
    // CIGAR operation Test
    if (Pv[bdp_idx] & mask) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (Mv[(bdp_idx-num_words64)] & mask) {
      operations[op_sentinel--] = 'I';
      --h;
    } else if ((text[h]==pattern[v])) {
      operations[op_sentinel--] = 'M';
      --h;
      --v;
    } else {
      operations[op_sentinel--] = 'X';
      --h;
      --v;
    }
  }
  while (h>=0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>=0) {operations[op_sentinel--] = 'D'; --v;}
  bpm_matrix->cigar->begin_offset = op_sentinel+1;
}
void bpm_compute(
    bpm_matrix_t* const bpm_matrix,
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const int text_length,
    const int max_distance) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix(
      bpm_matrix,bpm_pattern,
      text,text_length,max_distance);
  // Check distance
  if (bpm_matrix->min_score == UINT64_MAX) return;
  // Backtrace and generate CIGAR
  bpm_backtrace_matrix(bpm_matrix,bpm_pattern,text);
}
