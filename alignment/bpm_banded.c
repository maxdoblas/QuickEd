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
#include "utils/dna_text.h"

/*
 * Constants
 */
#define BPM_ALPHABET_LENGTH 4
#define BPM_W64_LENGTH UINT64_LENGTH
#define BPM_W64_SIZE   UINT64_SIZE
#define BPM_W64_ONES   UINT64_MAX
#define BPM_W64_MASK   (1ull<<63)

/*
 * Pattern Accessors
 */
#define banded_pattern_PEQ_IDX(word_pos,encoded_character)   ((word_pos*BPM_ALPHABET_LENGTH)+(encoded_character))
#define banded_pattern_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))
/*
 * Advance block functions (Improved)
 *   const @vector Eq,mask;
 *   return (Pv,Mv,PHout,MHout);
 */
#define BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout) \
  /* Computes modulator vector {Xv,Xh} ( cases A&C ) */ \
  const uint64_t Xv = Eq | Mv; \
  const uint64_t _Eq = Eq | MHin; \
  const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq; \
  /* Calculate Hout */ \
  uint64_t Ph = Mv | ~(Xh | Pv); \
  uint64_t Mh = Pv & Xh; \
  /* Account Hout that propagates for the next block */ \
  PHout = (Ph & mask)!=0; \
  MHout = (Mh & mask)!=0; \
  /* Hout become the Hin of the next cell */ \
  Ph <<= 1; \
  Mh <<= 1; \
  /* Account Hin coming from the previous block */ \
  Ph |= PHin; \
  Mh |= MHin; \
  /* Finally, generate the Vout */ \
  Pv = Mh | ~(Xv | Ph); \
  Mv = Ph & Xv

/*
 * Setup
 */
void banded_pattern_compile(
    banded_pattern_t* const banded_pattern,
    char* const pattern,
    const int pattern_length,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  const uint64_t PEQ_length = pattern_num_words64*BPM_W64_LENGTH;
  const uint64_t pattern_mod = pattern_length%BPM_W64_LENGTH;
  // Init fields
  banded_pattern->pattern = pattern;
  banded_pattern->pattern_length = pattern_length;
  banded_pattern->pattern_num_words64 = pattern_num_words64;
  banded_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words64*BPM_W64_SIZE;
  const uint64_t PEQ_size = BPM_ALPHABET_LENGTH*aux_vector_size;
  const uint64_t score_size = pattern_num_words64*UINT64_SIZE;
  const uint64_t total_memory = PEQ_size + 3*aux_vector_size + 2*score_size + (pattern_num_words64+1)*UINT64_SIZE;
  void* memory = mm_allocator_malloc(mm_allocator,total_memory);
  banded_pattern->PEQ = memory; memory += PEQ_size;
  banded_pattern->P = memory; memory += aux_vector_size;
  banded_pattern->M = memory; memory += aux_vector_size;
  banded_pattern->level_mask = memory; memory += aux_vector_size;
  banded_pattern->score = memory; memory += score_size;
  banded_pattern->init_score = memory; memory += score_size;
  banded_pattern->pattern_left = memory;
  // Init PEQ
  memset(banded_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = dna_encode(pattern[i]);
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    banded_pattern->PEQ[banded_pattern_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    uint64_t j;
    for (j=0;j<BPM_ALPHABET_LENGTH;++j) {
      banded_pattern->PEQ[banded_pattern_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words64-1;
  memset(banded_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    banded_pattern->level_mask[i] = BPM_W64_MASK;
    banded_pattern->init_score[i] = BPM_W64_LENGTH;
    banded_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  for (;i<=pattern_num_words64;++i) {
    banded_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  if (pattern_mod > 0) {
    const uint64_t mask_shift = pattern_mod-1;
    banded_pattern->level_mask[top] = 1ull<<(mask_shift);
    banded_pattern->init_score[top] = pattern_mod;
  } else {
    banded_pattern->level_mask[top] = BPM_W64_MASK;
    banded_pattern->init_score[top] = BPM_W64_LENGTH;
  }
}
void banded_pattern_free(
    banded_pattern_t* const banded_pattern,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_pattern->PEQ);
}
void banded_matrix_allocate(
    banded_matrix_t* const banded_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  banded_matrix->Mv = Mv;
  banded_matrix->Pv = Pv;
  // CIGAR
  banded_matrix->cigar = cigar_new(pattern_length+text_length);
  banded_matrix->cigar->end_offset = pattern_length+text_length;
}
void banded_matrix_free(
    banded_matrix_t* const banded_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_matrix->Mv);
  mm_allocator_free(mm_allocator,banded_matrix->Pv);
  // CIGAR
  cigar_free(banded_matrix->cigar);
}
/*
 * Edit distance computation using BPM
 */
void bpm_reset_search(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score) {
  // Reset score,P,M
  uint64_t i;
  P[0]=BPM_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    P[i]=BPM_W64_ONES;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }
}
void bpm_compute_matrix_banded(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length,
    const int bandwidth,
    const int max_distance) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const uint64_t num_words64 = banded_pattern->pattern_num_words64;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const level_mask = banded_pattern->level_mask;
  int64_t* const score = banded_pattern->score;
  const int64_t* const init_score = banded_pattern->init_score;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;
  bpm_reset_search(num_words64,Pv,Mv,score,init_score);

  const int k_end = ABS(text_length-pattern_length)+1;
  const int effective_bandwidth = MAX(k_end,bandwidth);

  // Advance in DP-bit_encoded matrix
  int lo_tmp = 0;
  int lo = 0;
  int hi_tmp = effective_bandwidth-1;
  int hi = 0;
  int last_hi = 0;
  uint64_t text_position;

  // Prologue: lo_band region
  for (text_position=0;text_position<=effective_bandwidth;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);

    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=lo;i<=last_hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = banded_pattern_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[banded_pattern_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    // Epilogue: Out of band adjacent blocks
    for (;i<=hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = banded_pattern_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = BPM_W64_ONES;
      uint64_t Mv_in = 0;
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[banded_pattern_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    last_hi = hi;
    // Compute new lo&hi limit
    hi_tmp++;
    hi = MIN(num_words64-1,hi_tmp/BPM_W64_SIZE);
  }

  // Main loop
  for (;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);

    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=lo;i<=last_hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = banded_pattern_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[banded_pattern_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    // Epilogue: Out of band adjacent blocks
    for (;i<=hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = banded_pattern_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = BPM_W64_ONES;
      uint64_t Mv_in = 0;
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[banded_pattern_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    last_hi = hi;
    // Compute lo&hi limit
    lo_tmp++;
    lo = lo_tmp/BPM_W64_LENGTH;
    hi_tmp++;
    hi = MIN(num_words64-1,hi_tmp/BPM_W64_SIZE);
  }

  // Return optimal column/distance
  const int64_t current_score = score[num_words64-1];
  if (current_score <= max_distance) {
    banded_matrix->min_score = current_score;
    banded_matrix->min_score_column = text_length-1;
  } else {
    banded_matrix->min_score = UINT64_MAX;
    banded_matrix->min_score_column = UINT64_MAX;
  }
}
void banded_backtrace_matrix(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = banded_pattern->pattern_num_words64;
  int64_t h = banded_matrix->min_score_column;
  int64_t v = pattern_length - 1;
  while (v >= 0 && h >= 0) {
    const uint8_t block = v / UINT64_LENGTH;
    const uint64_t bdp_idx = banded_pattern_BDP_IDX(h+1,num_words64,block);
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
  banded_matrix->cigar->begin_offset = op_sentinel+1;
}

void banded_compute(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length,
    const int bandwidth,
    const int max_distance) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded(
      banded_matrix,banded_pattern,
      text,text_length,bandwidth,max_distance);
  // Check distance
  if (banded_matrix->min_score == UINT64_MAX) return;
  // Backtrace and generate CIGAR
  banded_backtrace_matrix(banded_matrix,banded_pattern,text);
}
