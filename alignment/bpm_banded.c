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
#include <immintrin.h>

/*
 * Constants
 */
#define BPM_ALPHABET_LENGTH 5
#define BPM_W64_LENGTH UINT64_LENGTH
#define BPM_W64_SIZE   UINT64_SIZE
#define BPM_W64_ONES   UINT64_MAX
#define BPM_W64_MASK   (1ull<<63)

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   (((word_pos)*BPM_ALPHABET_LENGTH)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))
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

#define BPM_ADVANCE_BLOCK_LCS(Eq,mask,V,Hin,Hout) \
  const uint64_t _Eq = Eq; \
  const uint64_t H = V & Eq; \
  /* Account Hout that propagates for the next block */ \
  Hout = (H & mask)!=0; \
  /* Hout become the Hin of the next cell */ \
  H <<= 1; \
  /* Account Hin coming from the previous block */ \
  H |= Hin; \
  /* Finally, generate the Vout */ \
  V = Mh | ~(Xv | Ph); \

#define BPM_ADVANCE_BLOCK_NO_MASK(Eq,Pv,Mv,PHin,MHin,PHout,MHout) \
  /* Computes modulator vector {Xv,Xh} ( cases A&C ) */ \
  const uint64_t Xv = Eq | Mv; \
  const uint64_t _Eq = Eq | MHin; \
  const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq; \
  /* Calculate Hout */ \
  uint64_t Ph = Mv | ~(Xh | Pv); \
  uint64_t Mh = Pv & Xh; \
  /* Account Hout that propagates for the next block */ \
  PHout = Ph >> 63; \
  MHout = Mh >> 63; \
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
  const uint64_t total_memory = PEQ_size + 3*aux_vector_size + (pattern_num_words64)*UINT64_SIZE;
  void* memory = mm_allocator_malloc(mm_allocator,total_memory);
  banded_pattern->PEQ = memory; memory += PEQ_size;
  banded_pattern->P = memory; memory += aux_vector_size;
  banded_pattern->M = memory; memory += aux_vector_size;
  banded_pattern->level_mask = memory;
  // Init PEQ
  memset(banded_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = dna_encode(pattern[i]);
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    banded_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    uint64_t j;
    for (j=0;j<BPM_ALPHABET_LENGTH;++j) {
      banded_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  const uint64_t top = pattern_num_words64-1;
  memset(banded_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    banded_pattern->level_mask[i] = BPM_W64_MASK;
  }
  if (pattern_mod > 0) {
    const uint64_t mask_shift = pattern_mod-1;
    banded_pattern->level_mask[top] = 1ull<<(mask_shift);
  } else {
    banded_pattern->level_mask[top] = BPM_W64_MASK;
  }
}
void banded_pattern_free(
    banded_pattern_t* const banded_pattern,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_pattern->PEQ);
}

void banded_matrix_allocate_unaligned(
    banded_matrix_t* const banded_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    const int bandwidth,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const int k_end = ABS(text_length-pattern_length)+1;
  const int real_bandwidth = MAX(k_end,bandwidth);
  banded_matrix->effective_bandwidth_blocks = (real_bandwidth-1)/BPM_W64_LENGTH + 1;
  banded_matrix->effective_bandwidth = banded_matrix->effective_bandwidth_blocks*BPM_W64_LENGTH;

  const uint64_t num_words64 = banded_matrix->effective_bandwidth_blocks*2;
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
void banded_matrix_allocate_blocking(
    banded_matrix_t* const banded_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    const int bandwidth,
    mm_allocator_t* const mm_allocator){
  // Parameters
  const int k_end = ABS(text_length-pattern_length)+1;
  const int real_bandwidth = MAX(k_end,bandwidth);
  banded_matrix->effective_bandwidth_blocks = DIV_CEIL(real_bandwidth, BPM_W64_LENGTH) + 1;
  banded_matrix->effective_bandwidth = banded_matrix->effective_bandwidth_blocks*BPM_W64_LENGTH;

  const uint64_t num_words64 = banded_matrix->effective_bandwidth_blocks;
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

void banded_matrix_allocate_cutoff(
    banded_matrix_t* const banded_matrix,
    const int64_t pattern_length,
    const int64_t text_length,
    const int64_t cutoff_score,
    mm_allocator_t* const mm_allocator){

  const int64_t k_end = ABS(((int64_t)text_length)-(int64_t)(pattern_length))+1;
  banded_matrix->cutoff_score = MAX(MAX(k_end,cutoff_score),65);
  banded_matrix->sequence_length_diff = pattern_length - text_length;
  banded_matrix->relative_cutoff_score = DIV_CEIL((banded_matrix->cutoff_score - ABS(banded_matrix->sequence_length_diff)), 2);
  if(banded_matrix->sequence_length_diff >= 0){
    banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH);
    banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score + banded_matrix->sequence_length_diff, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
  } else {
    banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score - banded_matrix->sequence_length_diff, BPM_W64_LENGTH);
    banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
  }
  banded_matrix->effective_bandwidth = banded_matrix->cutoff_score;

  const int64_t num_words64 = banded_matrix->effective_bandwidth_blocks;
  // Allocate auxiliary matrix
  const int64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const scores = (uint64_t*)mm_allocator_malloc(mm_allocator, (DIV_CEIL(pattern_length, BPM_W64_LENGTH) + num_words64/2) * UINT64_SIZE);
  banded_matrix->Mv = Mv;
  banded_matrix->Pv = Pv;
  banded_matrix->scores = scores;
  // CIGAR
  banded_matrix->cigar = cigar_new(pattern_length+text_length);
  banded_matrix->cigar->end_offset = pattern_length+text_length;
}

void banded_matrix_allocate_cutoff_score(
    banded_matrix_t* const banded_matrix,
    const int64_t pattern_length,
    const int64_t text_length,
    const int64_t cutoff_score,
    mm_allocator_t* const mm_allocator){
  // Parameters
  const int64_t k_end = ABS(((int64_t)text_length)-(int64_t)(pattern_length))+1;
  banded_matrix->cutoff_score = MAX(MAX(k_end,cutoff_score),65);
  banded_matrix->sequence_length_diff = pattern_length - text_length;
  banded_matrix->relative_cutoff_score = DIV_CEIL((banded_matrix->cutoff_score - ABS(banded_matrix->sequence_length_diff)), 2);
  if(banded_matrix->sequence_length_diff >= 0){
    banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH);
    banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score + banded_matrix->sequence_length_diff, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
  } else {
    banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score - banded_matrix->sequence_length_diff, BPM_W64_LENGTH);
    banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
  }
  banded_matrix->effective_bandwidth = banded_matrix->cutoff_score;

  const uint64_t num_words64 = banded_matrix->effective_bandwidth_blocks;
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE; /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const scores = (uint64_t*)mm_allocator_malloc(mm_allocator, (DIV_CEIL(pattern_length, BPM_W64_LENGTH) + num_words64/2) * UINT64_SIZE);
  banded_matrix->Mv = Mv;
  banded_matrix->Pv = Pv;
  banded_matrix->scores = scores;
  // CIGAR
  banded_matrix->cigar = cigar_new(pattern_length+text_length);
  banded_matrix->cigar->end_offset = pattern_length+text_length;
}

void banded_matrix_allocate(
    banded_matrix_t* const banded_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    const int bandwidth,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  banded_matrix->Mv = Mv;
  banded_matrix->Pv = Pv;

  const int k_end = ABS(text_length-pattern_length)+1;
  const int effective_bandwidth = MAX(k_end,bandwidth);
  banded_matrix->effective_bandwidth = effective_bandwidth;
  // Lower and upper bounds
  int* const lo = (int*)mm_allocator_malloc(mm_allocator,(text_length+1)*sizeof(int));
  int* const hi = (int*)mm_allocator_malloc(mm_allocator,(text_length+1)*sizeof(int));
  banded_matrix->lo = lo;
  banded_matrix->hi = hi;

  uint64_t lo_tmp = 0;
  uint64_t hi_tmp = effective_bandwidth-1;

  int i;
  for (i = 0; i <= effective_bandwidth; i++) {
    lo[i] = 0;
    hi_tmp++;
    hi[i] = MIN(num_words64-1,hi_tmp/BPM_W64_LENGTH);
  }
  for (; (hi[i-1]) < (num_words64-1); i++) {
    lo_tmp++;
    lo[i] = lo_tmp/BPM_W64_LENGTH;
    hi_tmp++;
    hi[i] = MIN(num_words64-1,hi_tmp/BPM_W64_LENGTH);
  }
  for (; i <= text_length; i++) {
    lo_tmp++;
    lo[i] = lo_tmp/BPM_W64_LENGTH;
    hi[i] = num_words64-1;
  }

  // CIGAR
  banded_matrix->cigar = cigar_new(pattern_length+text_length);
  banded_matrix->cigar->end_offset = pattern_length+text_length;
}
void banded_matrix_free(
    banded_matrix_t* const banded_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_matrix->Mv);
  mm_allocator_free(mm_allocator,banded_matrix->Pv);
  mm_allocator_free(mm_allocator,banded_matrix->lo);
  mm_allocator_free(mm_allocator,banded_matrix->hi);
  // CIGAR
  cigar_free(banded_matrix->cigar);
}

void banded_matrix_free_unaligned(
    banded_matrix_t* const banded_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_matrix->Mv);
  mm_allocator_free(mm_allocator,banded_matrix->Pv);
  // CIGAR
  cigar_free(banded_matrix->cigar);
}

void banded_matrix_free_blocking(
    banded_matrix_t* const banded_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_matrix->Mv);
  mm_allocator_free(mm_allocator,banded_matrix->Pv);
  // CIGAR
  cigar_free(banded_matrix->cigar);
}

void banded_matrix_free_cutoff(
    banded_matrix_t* const banded_matrix,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,banded_matrix->Mv);
  mm_allocator_free(mm_allocator,banded_matrix->Pv);
  mm_allocator_free(mm_allocator,banded_matrix->scores);
  // CIGAR
  cigar_free(banded_matrix->cigar);
}
/*
 * Edit distance computation using BPM
 */
void bpm_reset_search(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M) {
  // Reset P,M
  uint64_t i;
  P[0]=BPM_W64_ONES;
  M[0]=0;
  for (i=1;i<num_words;++i) {
    P[i]=BPM_W64_ONES;
    M[i]=0;
  }
}

void bpm_reset_search_banded_score(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M,
    uint64_t* const scores) {
  // Reset P,M
  uint64_t i;
  P[0]=BPM_W64_ONES;
  M[0]=0;
  scores[0] = BPM_W64_LENGTH;
  for (i=1;i<num_words;++i) {
    P[i]=BPM_W64_ONES;
    M[i]=0;
    scores[i] = scores[i-1] + BPM_W64_LENGTH;
  }
}

void bpm_compute_matrix_banded(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const uint64_t num_words64 = banded_pattern->pattern_num_words64;
  //const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;
  //const int effective_bandwidth = banded_matrix->effective_bandwidth;
  const int* lo_vec = banded_matrix->lo;
  const int* hi_vec = banded_matrix->hi;
  bpm_reset_search(num_words64,Pv,Mv);

  // Advance in DP-bit_encoded matrix
  int last_hi = 0;
  uint64_t text_position;
  //uint64_t count = 0;

  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);

    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;

    // Main Loop
    const int lo = lo_vec[text_position];
    const int hi = hi_vec[text_position];
    for (i=lo;i<=last_hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    // Epilogue: Out of band adjacent blocks
    for (;i<=hi;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = BPM_W64_ONES;
      uint64_t Mv_in = 0;
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }

    last_hi = hi;
  }
  //printf("Number of BPM_BLOCKS = %d\n",count);

}
void banded_backtrace_matrix(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = banded_pattern->pattern_num_words64;
  int64_t h = text_length - 1;
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
  banded_matrix->cigar->begin_offset = op_sentinel+1;
}

void bpm_compute_matrix_banded_unaligned(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  const int effective_bandwidth = banded_matrix->effective_bandwidth;
  const uint64_t num_words64 = effective_bandwidth_blocks*2;
  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;

  bpm_reset_search(effective_bandwidth_blocks*2,Pv,Mv);

  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  //uint64_t count = 0;
  uint64_t prologue_rows = MIN(text_length,(effective_bandwidth-1));
  // Prologue: lo_band region
  for (text_position=0;text_position<prologue_rows;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=0;i<effective_bandwidth_blocks*2;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }
  }
  uint64_t pos_v = 0;
  // Main loop
  for (;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    const uint64_t shift = pos_v % BPM_W64_LENGTH;
    const uint64_t pos_v_block = pos_v / BPM_W64_LENGTH;
    const uint64_t shift_mask = shift ? 0xFFFFFFFFFFFFFFFFULL : 0ULL;

    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=0;i<effective_bandwidth_blocks*2;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v_block,enc_char)] >> shift | ((PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v_block+1,enc_char)] << (BPM_W64_LENGTH - shift)) & shift_mask);

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;
      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in >> 1;
      Mv[next_bdp_idx] = Mv_in >> 1;
      Pv[next_bdp_idx - 1] |= (uint64_t)((Pv_in & 0x1UL) || i==0) << 63;
      Mv[next_bdp_idx - 1] |= (uint64_t)((Mv_in & 0x1UL) && i!=0) << 63;
      PHin=PHout;
      MHin=MHout;
    }
    pos_v++;
  }
}


void banded_backtrace_matrix_unaligned(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  const int effective_bandwidth = banded_matrix->effective_bandwidth;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = 2*effective_bandwidth_blocks;
  int64_t h = text_length - 1;
  int64_t v = pattern_length - 1;

  while (v >= 0 && h >= 0) {
    const uint64_t effective_v = (h >= effective_bandwidth - 1) ? (v - (h + 1 - effective_bandwidth)) : v;
    const uint64_t effective_v_r = (h >= effective_bandwidth - 1) ? (v - (h + 1 - effective_bandwidth) - 1) : v;
    const uint64_t block = effective_v / UINT64_LENGTH;
    const uint64_t block_r = effective_v_r / UINT64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h,num_words64,block);
    const uint64_t bdp_idx_r = BPM_PATTERN_BDP_IDX(h+1,num_words64,block_r);
    const uint64_t mask = 1UL << (effective_v % UINT64_LENGTH);
    const uint64_t mask_r = 1UL << (effective_v_r % UINT64_LENGTH);

    // CIGAR operation Test
    if (Pv[bdp_idx_r] & mask_r) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (Mv[(bdp_idx)] & mask) {
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

void bpm_compute_matrix_banded_blocking(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  const uint64_t num_words64 = effective_bandwidth_blocks;
  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;

  bpm_reset_search(effective_bandwidth_blocks,Pv,Mv);

  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  //uint64_t count = 0;

  uint64_t prologue_rows = MIN(text_length,(effective_bandwidth_blocks/2)*BPM_W64_LENGTH);
  // Prologue: lo_band region
  for (text_position=0;text_position<prologue_rows;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=0;i<effective_bandwidth_blocks-1;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }
  }

  Pv[BPM_PATTERN_BDP_IDX(text_position,num_words64,effective_bandwidth_blocks-1)] = BPM_W64_ONES;
  Mv[BPM_PATTERN_BDP_IDX(text_position,num_words64,effective_bandwidth_blocks-1)] = 0;

  uint64_t pos_v = 0;
  // Main loop
  for (;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i,PHin=1,MHin=0,PHout,MHout;
    // Main Loop
    for (i=0;i<effective_bandwidth_blocks;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i+pos_v];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
    }
    
    // Shift results one block in the last column of a 64-column block
    if((text_position + 1) % 64 == 0){
      uint64_t next_bdp_idx = BPM_PATTERN_BDP_IDX(text_position+1,num_words64,0);
      for(int j = 0; j < effective_bandwidth_blocks-1; j++){
        Pv[next_bdp_idx] = Pv[next_bdp_idx+1];
        Mv[next_bdp_idx] = Mv[next_bdp_idx+1];
        next_bdp_idx++;
      }
      Pv[next_bdp_idx] = BPM_W64_ONES;
      Mv[next_bdp_idx] = 0;
      pos_v++;
    } 
  }
  //printf("Number of BPM_BLOCKS = %d\n",count);
}


void banded_backtrace_matrix_blocking(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = effective_bandwidth_blocks;
  int64_t h = text_length - 1;
  int64_t v = pattern_length - 1;

  const uint64_t half_band_blocks = effective_bandwidth_blocks / 2;
  while (v >= 0 && h >= 0) {
    const uint64_t block_h = h / BPM_W64_LENGTH;
    const uint64_t block_h_r = (h+1) / BPM_W64_LENGTH;
    const uint64_t effective_v = (block_h > half_band_blocks) ? v - BPM_W64_LENGTH*(block_h-half_band_blocks) : v;
    const uint64_t effective_v_r = (block_h_r > half_band_blocks) ? v - BPM_W64_LENGTH*(block_h_r-half_band_blocks) : v;
    const uint64_t block_v = effective_v / BPM_W64_LENGTH;
    const uint64_t block_v_r = effective_v_r / BPM_W64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h,num_words64,block_v);
    const uint64_t bdp_idx_r = BPM_PATTERN_BDP_IDX(h+1,num_words64,block_v_r);
    const uint64_t mask = 1UL << (effective_v % BPM_W64_LENGTH);
    const uint64_t mask_r = 1UL << (effective_v_r % BPM_W64_LENGTH);

    // CIGAR operation Test
    if (Pv[bdp_idx_r] & mask_r) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (Mv[(bdp_idx)] & mask) {
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

void bpm_compute_matrix_banded_cutoff(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int64_t text_length,
    const int64_t cutoff_score) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const int64_t effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;

  const int64_t num_block_rows = DIV_CEIL(banded_pattern->pattern_length,BPM_W64_LENGTH);

  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;
  uint64_t* const scores = banded_matrix->scores;
  const uint64_t num_words64 = effective_bandwidth_blocks;

  const int64_t sequence_length_diff = banded_matrix->sequence_length_diff;
  const int64_t prologue_columns = banded_matrix->prolog_column_blocks;
  const int64_t finish_v_pos_inside_band = prologue_columns * BPM_W64_LENGTH + sequence_length_diff;

  // Prepare last block of the next column
  int64_t pos_v = -prologue_columns;
  int64_t pos_h = 0;
  int64_t first_block_v = prologue_columns;
  int64_t last_block_v = effective_bandwidth_blocks-1;

  bpm_reset_search_banded_score(effective_bandwidth_blocks,Pv,Mv,scores);

  // Advance in DP-bit_encoded matrix
  int64_t text_position;
  // Main loop
  for (text_position = 0; text_position < text_length; ++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    int64_t i;
    uint64_t PHin = 1, MHin = 0, PHout, MHout;
    // Main Loop
    for (i = first_block_v; i <= last_block_v; ++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i+pos_v];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX((i+pos_v),enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;

      PHin=PHout;
      MHin=MHout;
      scores[i+pos_v] = scores[i+pos_v] + PHout - MHout;
    }
    
    // Update the score for the new column
    if((text_position + 1) % 64 == 0){
      //printf("-----------------------------------------------------\n");
      // chech if the band of the lower side should be cutted
      int cut_band_lower = (first_block_v+2 < last_block_v) && (finish_v_pos_inside_band > BPM_W64_LENGTH*(first_block_v+1))
                           && (scores[first_block_v + pos_v + 1] + (finish_v_pos_inside_band - BPM_W64_LENGTH*(first_block_v+1))) > banded_matrix->cutoff_score;

      if (cut_band_lower && (pos_h >= prologue_columns)){
          first_block_v++;
      }else if (!cut_band_lower && (pos_h < prologue_columns)){
          first_block_v--;
      }
      
      uint64_t next_bdp_idx = BPM_PATTERN_BDP_IDX(text_position+1,num_words64,0);
      // Shift results one block in the last column of a 64-column block
      for(int64_t j = first_block_v; j < last_block_v; j++){
        Pv[next_bdp_idx+j] = Pv[next_bdp_idx+j+1];
        Mv[next_bdp_idx+j] = Mv[next_bdp_idx+j+1];
      }
      Pv[next_bdp_idx+last_block_v] = BPM_W64_ONES;
      Mv[next_bdp_idx+last_block_v] = 0;
      // Update the score for the new column
      int64_t pos = last_block_v + pos_v;
      scores[pos+1] = scores[pos] + BPM_W64_LENGTH;

      // chech if the band of the higher side should be cutted
      int cut_band_higer = (first_block_v+2 < last_block_v) && (BPM_W64_LENGTH*(last_block_v-1) > finish_v_pos_inside_band)
                           && (scores[last_block_v + pos_v - 1] + (BPM_W64_LENGTH*(last_block_v-1)- finish_v_pos_inside_band)) > banded_matrix->cutoff_score;

      if (cut_band_higer || (pos_v+last_block_v >= num_block_rows-1)){
        last_block_v--;
      }
      pos_v++; pos_h++;  
    } 
  }

  uint64_t final_score;
  if (banded_pattern->pattern_length%BPM_W64_LENGTH){
    final_score = scores[(banded_pattern->pattern_length)/BPM_W64_LENGTH]-(BPM_W64_LENGTH-(banded_pattern->pattern_length%BPM_W64_LENGTH));
  }else{
    final_score = scores[(banded_pattern->pattern_length-1)/BPM_W64_LENGTH];
  }
  banded_matrix->cigar->score = final_score;
  banded_matrix->higher_block = last_block_v;
  banded_matrix->lower_block = first_block_v;
}

void bpm_compute_matrix_banded_cutoff_score(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int64_t text_length,
    const int64_t text_finish_pos,
    const int64_t cutoff_score) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  // TODO: remove if necessary
  const int64_t k_end = ABS(((int64_t)text_length)-(int64_t)(banded_pattern->pattern_length))+1;
  const int64_t real_bandwidth = MAX(MAX(k_end,cutoff_score),65);
  const int64_t effective_bandwidth_blocks = DIV_CEIL(real_bandwidth, BPM_W64_LENGTH) + 1;
  const int64_t num_block_rows = DIV_CEIL(banded_pattern->pattern_length,BPM_W64_LENGTH);

  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;
  uint64_t* const scores = banded_matrix->scores;

  const int64_t sequence_length_diff = banded_matrix->sequence_length_diff;
  const int64_t prologue_columns = banded_matrix->prolog_column_blocks;
  const int64_t finish_v_pos_inside_band = prologue_columns * BPM_W64_LENGTH + sequence_length_diff;

  // Prepare last block of the next column
  int64_t pos_v = -prologue_columns;
  int64_t pos_h = 0;
  int64_t first_block_v = prologue_columns;
  int64_t last_block_v = effective_bandwidth_blocks-1;

  bpm_reset_search_banded_score(effective_bandwidth_blocks,Pv,Mv,scores);

  // Advance in DP-bit_encoded matrix
  int64_t text_position;
  //uint64_t count = 0;
  // Main loop
  for (text_position = 0; text_position < text_finish_pos; ++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    int64_t i;
    uint64_t PHin = 1, MHin = 0, PHout, MHout;
    // Main Loop
    for (i = first_block_v; i <= last_block_v; ++i) {
      /* Calculate Step Data */
      uint64_t Pv_in = Pv[i];
      uint64_t Mv_in = Mv[i];
      const uint64_t mask = level_mask[i+pos_v];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX((i+pos_v),enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);

      /* Swap propagate Hv */
      Pv[i] = Pv_in;
      Mv[i] = Mv_in;
      PHin=PHout;
      MHin=MHout;
      scores[i+pos_v] = scores[i+pos_v] + PHout - MHout;
    }
    
    // Update the score for the new column
    
    if((text_position + 1) % 64 == 0){

      // chech if the band of the lower side should be cutted
      int cut_band_lower = (first_block_v+2 < last_block_v) && (finish_v_pos_inside_band > BPM_W64_LENGTH*(first_block_v+1))
                           && (scores[first_block_v + pos_v + 1] + (finish_v_pos_inside_band - BPM_W64_LENGTH*(first_block_v+1))) > banded_matrix->cutoff_score;

      // if we are in the prolog columns, we have to decrease the first_block_v (apart of cutting the band if necessary)
      if (cut_band_lower && (pos_h >= prologue_columns)){
          first_block_v++;
      }else if (!cut_band_lower && (pos_h < prologue_columns)){
          first_block_v--;
      }

      // Shift results one block in the last column of a 64-column block
      for(int64_t j = first_block_v; j < last_block_v; j++){
        Pv[j] = Pv[j+1];
        Mv[j] = Mv[j+1];
      }
      Pv[last_block_v] = BPM_W64_ONES;
      Mv[last_block_v] = 0;
      // Update the score for the new column
      int64_t pos = last_block_v + pos_v;
      scores[pos+1] = scores[pos] + BPM_W64_LENGTH;

      // chech if the band of the higher side should be cutted
      int cut_band_higer = (first_block_v+2 < last_block_v) && (BPM_W64_LENGTH*(last_block_v-1) > finish_v_pos_inside_band)
                           && (scores[last_block_v + pos_v - 1] + (BPM_W64_LENGTH*(last_block_v-1)- finish_v_pos_inside_band)) > banded_matrix->cutoff_score;

      if (cut_band_higer || (pos_v+last_block_v >= num_block_rows)){
        last_block_v--;
      }
      pos_v++; pos_h++;     
    } 
  }
  uint64_t final_score;
  if (banded_pattern->pattern_length%BPM_W64_LENGTH){
    final_score = scores[(banded_pattern->pattern_length)/BPM_W64_LENGTH]-(BPM_W64_LENGTH-(banded_pattern->pattern_length%BPM_W64_LENGTH));
  }else{
    final_score = scores[(banded_pattern->pattern_length-1)/BPM_W64_LENGTH];
  }
  banded_matrix->cigar->score = final_score;
  banded_matrix->higher_block = last_block_v;
  banded_matrix->lower_block = first_block_v;
}

void bpm_compute_matrix_banded_cutoff_simd(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length,
    const uint64_t cutoff_score) {
  // Pattern variables
  const uint64_t* PEQ = banded_pattern->PEQ;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  const uint64_t num_words64 = effective_bandwidth_blocks;
  const uint64_t* const level_mask = banded_pattern->level_mask;
  uint64_t* const Pv = banded_matrix->Pv;
  uint64_t* const Mv = banded_matrix->Mv;
  uint64_t* const scores = banded_matrix->scores;

  const uint64_t finish_v_pos_inside_band = (effective_bandwidth_blocks/2) * BPM_W64_LENGTH - (text_length - banded_pattern->pattern_length);
  bpm_reset_search_banded_score(effective_bandwidth_blocks,Pv,Mv,scores);

  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  //uint64_t count = 0;
  uint64_t prologue_rows = MIN(text_length,(effective_bandwidth_blocks/2)*BPM_W64_LENGTH);
  // Prologue: lo_band region
  for (text_position = 0; text_position < prologue_rows; ++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i, PHin = 1, MHin = 0, PHout, MHout;
    // Main Loop
    for (i = 0; i < (effective_bandwidth_blocks - 1); ++i) {

      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;

      PHin=PHout;
      MHin=MHout;
      scores[i] = scores[i] + PHout - MHout;
    }
  }

  // Prepare last block of the next column
  uint64_t pos = effective_bandwidth_blocks-1;
  scores[pos] = scores[pos-1] + BPM_W64_LENGTH;
  Pv[BPM_PATTERN_BDP_IDX(text_position,num_words64,effective_bandwidth_blocks-1)] = BPM_W64_ONES;
  Mv[BPM_PATTERN_BDP_IDX(text_position,num_words64,effective_bandwidth_blocks-1)] = 0;

  uint64_t pos_v = 0;
  uint64_t first_block_v = 0;
  uint64_t last_block_v = effective_bandwidth_blocks-1;
  // Main loop
  //uint64_t last_text = text_length - effective_bandwidth;
  for (;text_position<text_length-63;text_position+=64) {
    // Advance all blocks
    // Advance all blocks
    uint64_t i, PHin[65], MHin[65];
    for (uint64_t text_pos_i=0;text_pos_i<65;text_pos_i++){
      PHin[text_pos_i] = 1;
      MHin[text_pos_i] = 0;
    }
    // Main Loop SSE (antidiagonal parallelisem)
    //  1 2 3 4 5 ...  64
    //  2 3 4 5 ... 64 65
    __m128i PHout_v,MHout_v,PHin_v,MHin_v;
    for (i = first_block_v; i <= last_block_v-1; i+=2) {
      if ((level_mask[i+pos_v] | level_mask[i+pos_v+1]) != 0x8000000000000000ULL){
        break;
      }

      uint64_t Ph_first=PHin[0],Mh_first=MHin[0];
      {
        const uint8_t enc_char = dna_encode(text[text_position]);
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
        const uint64_t next_bdp_idx = bdp_idx+num_words64;
        uint64_t Pv_in = Pv[bdp_idx];
        uint64_t Mv_in = Mv[bdp_idx];
        uint64_t PHout,MHout;
        const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v,enc_char)];
        BPM_ADVANCE_BLOCK_NO_MASK(Eq,Pv_in,Mv_in,Ph_first,Mh_first,PHout,MHout);
        Pv[next_bdp_idx] = Pv_in;
        Mv[next_bdp_idx] = Mv_in;
        Ph_first = PHout;
        Mh_first = MHout;
        scores[i+pos_v] = scores[i+pos_v] + PHout - MHout;
      }
      PHin_v = _mm_set_epi64x(Ph_first,PHin[1]);
      MHin_v = _mm_set_epi64x(Mh_first,MHin[1]); 
      
      __m128i scores_v = _mm_set_epi64x(scores[i+pos_v+1],scores[i+pos_v]);

      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position+1,num_words64,i);
      const uint64_t bdp_idx_2 = BPM_PATTERN_BDP_IDX(text_position,num_words64,i+1);
      __m128i Mv_v = _mm_set_epi64x(Mv[bdp_idx_2], Mv[bdp_idx]);
      __m128i Pv_v = _mm_set_epi64x(Pv[bdp_idx_2], Pv[bdp_idx]);
      uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t next_bdp_idx_2 = bdp_idx_2+num_words64;

      for (uint64_t text_pos_i=1;text_pos_i<64;text_pos_i++) {
        // Fetch next character
        uint8_t enc_char = dna_encode(text[text_position+text_pos_i]);
        uint8_t enc_char_2 = dna_encode(text[text_position+text_pos_i-1]);
        
        /* Calculate Step Data */
        const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v,enc_char)];
        const uint64_t Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v+1,enc_char_2)];
        __m128i Eq_v = _mm_set_epi64x(Eq_2,Eq); 
        
        __m128i Xv =  _mm_or_si128(Eq_v, Mv_v); 
        __m128i _Eq = _mm_or_si128(Eq_v, MHin_v); 
        __m128i sum_128 = _mm_add_epi64(_mm_and_si128(_Eq, Pv_v), Pv_v);
        __m128i Xh =  _mm_or_si128(_mm_xor_si128(sum_128, Pv_v), _Eq); 
        /* Calculate Hout */ 
        __m128i Ph_v = _mm_or_si128(Mv_v, ~(_mm_or_si128(Xh, Pv_v))); 
        __m128i Mh_v = _mm_and_si128(Pv_v, Xh); 
        /* Account Hout that propagates for the next block */ 
        PHout_v = _mm_srli_epi64(Ph_v, 63); 
        MHout_v = _mm_srli_epi64(Mh_v, 63); 
        /* Hout become the Hin of the next cell */ 
        Ph_v = _mm_slli_epi64(Ph_v, 1);  
        Mh_v = _mm_slli_epi64(Mh_v, 1); 
        /* Account Hin coming from the previous block */ 
        Ph_v = _mm_or_si128(Ph_v,PHin_v); 
        Mh_v = _mm_or_si128(Mh_v,MHin_v); 
        /* Finally, generate the Vout */ 
        Pv_v = _mm_or_si128(Mh_v, ~(_mm_or_si128(Xv, Ph_v))); 
        Mv_v = _mm_and_si128(Ph_v, Xv);

        Pv[next_bdp_idx] = _mm_extract_epi64(Pv_v,0);
        Mv[next_bdp_idx] = _mm_extract_epi64(Mv_v,0);
        Pv[next_bdp_idx_2] = _mm_extract_epi64(Pv_v,1);
        Mv[next_bdp_idx_2] = _mm_extract_epi64(Mv_v,1);

        next_bdp_idx += num_words64;
        next_bdp_idx_2 += num_words64;

        scores_v = _mm_add_epi64(_mm_sub_epi64(scores_v,MHout_v),PHout_v);

        const uint64_t PHout_0 = _mm_extract_epi64(PHout_v,0);
        const uint64_t MHout_0 = _mm_extract_epi64(MHout_v,0);
        const uint64_t PHout_1 = _mm_extract_epi64(PHout_v,1);
        const uint64_t MHout_1 = _mm_extract_epi64(MHout_v,1);
        PHin_v=_mm_set_epi64x(PHout_0, PHin[text_pos_i+1]);
        MHin_v=_mm_set_epi64x(MHout_0, MHin[text_pos_i+1]);
        PHin[text_pos_i-1]=PHout_1;
        MHin[text_pos_i-1]=MHout_1;
      }
      scores[i+pos_v]=_mm_extract_epi64(scores_v,0);
      scores[i+pos_v+1]=_mm_extract_epi64(scores_v,1);
      
      Ph_first=_mm_extract_epi64(PHout_v,0);
      Mh_first=_mm_extract_epi64(MHout_v,0);
      {
        const uint8_t enc_char = dna_encode(text[text_position+63]);
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position+63,num_words64,i+1);
        const uint64_t next_bdp_idx = bdp_idx+num_words64;
        uint64_t Pv_in = Pv[bdp_idx];
        uint64_t Mv_in = Mv[bdp_idx];
        uint64_t PHout,MHout;
        const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v+1,enc_char)];
        BPM_ADVANCE_BLOCK_NO_MASK(Eq,Pv_in,Mv_in,Ph_first,Mh_first,PHout,MHout);
        Pv[next_bdp_idx] = Pv_in;
        Mv[next_bdp_idx] = Mv_in;

        Ph_first = PHout;
        Mh_first = MHout;
        scores[i+pos_v+1] = scores[i+pos_v+1] + PHout - MHout;
        PHin[63]=PHout;
        MHin[63]=MHout;
      }
    }

    for (; i <= last_block_v; ++i) {
      for (uint64_t text_pos_i=0;text_pos_i<64;text_pos_i++) {
        // Fetch next character
        const uint8_t enc_char = dna_encode(text[text_position+text_pos_i]);
        /* Calculate Step Data */
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position+text_pos_i,num_words64,i);
        const uint64_t next_bdp_idx = bdp_idx+num_words64;
        uint64_t Pv_in = Pv[bdp_idx];
        uint64_t Mv_in = Mv[bdp_idx];
        const uint64_t mask = level_mask[i+pos_v];
        const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v,enc_char)];
        uint64_t PHout, MHout;
        /* Compute Block */
        BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin[text_pos_i],MHin[text_pos_i],PHout,MHout);
        //count++;

        /* Swap propagate Hv */
        Pv[next_bdp_idx] = Pv_in;
        Mv[next_bdp_idx] = Mv_in;
        PHin[text_pos_i]=PHout;
        MHin[text_pos_i]=MHout;
        scores[i+pos_v] = scores[i+pos_v] + PHout - MHout;
      }
    }

    // Update the score for the new column
    // Shift results one block in the last column of a 64-column block
    uint64_t next_bdp_idx = BPM_PATTERN_BDP_IDX(text_position+64,num_words64,0);
    for(int j = first_block_v; j < last_block_v; j++){
      Pv[next_bdp_idx+j] = Pv[next_bdp_idx+j+1];
      Mv[next_bdp_idx+j] = Mv[next_bdp_idx+j+1];
    }
    Pv[next_bdp_idx+last_block_v] = BPM_W64_ONES;
    Mv[next_bdp_idx+last_block_v] = 0;
    
    uint64_t pos = last_block_v + pos_v;
    scores[pos+1] = scores[pos] + BPM_W64_LENGTH;

    if ((scores[first_block_v + pos_v + 1] + (finish_v_pos_inside_band - BPM_W64_LENGTH*(first_block_v + 1))) > cutoff_score && (first_block_v+2 < last_block_v) ){
      first_block_v++;
    }
    if ((scores[last_block_v + pos_v - 1] + (BPM_W64_LENGTH*(last_block_v - 1) - finish_v_pos_inside_band))> cutoff_score && (first_block_v+2 < last_block_v) ){
      last_block_v--;
    }
    pos_v++;
  }

  for (;text_position<text_length;text_position++) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i, PHin = 1, MHin = 0, PHout, MHout;
    // Main Loop
    for (i = first_block_v; i <= last_block_v; ++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i+pos_v];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i+pos_v,enc_char)];

      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      //count++;

      /* Swap propagate Hv */
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
      PHin=PHout;
      MHin=MHout;
      scores[i+pos_v] = scores[i+pos_v] + PHout - MHout;
    }
  }
  //printf("Number of BPM_BLOCKS = %d\n",count);
}

void banded_backtrace_matrix_cutoff(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text,
    const int64_t text_length,
    const int64_t cutoff_score) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  const int64_t prologue_columns = banded_matrix->prolog_column_blocks;

  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = effective_bandwidth_blocks;
  int64_t h = text_length - 1;
  int64_t v = pattern_length - 1;

  while (v >= 0 && h >= 0) {
    const int64_t block_h = h / BPM_W64_LENGTH;
    const int64_t block_h_r = (h+1) / BPM_W64_LENGTH;
    const int64_t effective_v = v - BPM_W64_LENGTH*(block_h - prologue_columns);
    const int64_t effective_v_r = v - BPM_W64_LENGTH*(block_h_r - prologue_columns);
    const int64_t block_v = effective_v / BPM_W64_LENGTH;
    const int64_t block_v_r = effective_v_r / BPM_W64_LENGTH;
    const int64_t bdp_idx = BPM_PATTERN_BDP_IDX(h,num_words64,block_v);
    const int64_t bdp_idx_r = BPM_PATTERN_BDP_IDX(h+1,num_words64,block_v_r);
    const uint64_t mask = 1UL << (effective_v % BPM_W64_LENGTH);
    const uint64_t mask_r = 1UL << (effective_v_r % BPM_W64_LENGTH);

    // CIGAR operation Test
    if (Pv[bdp_idx_r] & mask_r) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (Mv[(bdp_idx)] & mask) {
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

void banded_backtrace_matrix_cutoff_simd(
    banded_matrix_t* const banded_matrix,
    const banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Parameters
  char* const pattern = banded_pattern->pattern;
  const uint64_t pattern_length = banded_pattern->pattern_length;
  const uint64_t* const Pv = banded_matrix->Pv;
  const uint64_t* const Mv = banded_matrix->Mv;
  char* const operations = banded_matrix->cigar->operations;
  int op_sentinel = banded_matrix->cigar->end_offset-1;
  const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = effective_bandwidth_blocks;
  int64_t h = text_length - 1;
  int64_t v = pattern_length - 1;

  const uint64_t half_band_blocks = effective_bandwidth_blocks / 2;
  while (v >= 0 && h >= 0) {
    const uint64_t block_h = h / BPM_W64_LENGTH;
    const uint64_t block_h_r = (h+1) / BPM_W64_LENGTH;
    const uint64_t effective_v = (block_h > half_band_blocks) ? v - BPM_W64_LENGTH*(block_h-half_band_blocks) : v;
    const uint64_t effective_v_r = (block_h_r > half_band_blocks) ? v - BPM_W64_LENGTH*(block_h_r-half_band_blocks) : v;
    const uint64_t block_v = effective_v / BPM_W64_LENGTH;
    const uint64_t block_v_r = effective_v_r / BPM_W64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h,num_words64,block_v);
    const uint64_t bdp_idx_r = BPM_PATTERN_BDP_IDX(h+1,num_words64,block_v_r);
    const uint64_t mask = 1UL << (effective_v % BPM_W64_LENGTH);
    const uint64_t mask_r = 1UL << (effective_v_r % BPM_W64_LENGTH);
    // CIGAR operation Test
    if (Pv[bdp_idx_r] & mask_r) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (Mv[(bdp_idx)] & mask) {
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
    const int text_length) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded(
      banded_matrix,banded_pattern,
      text,text_length);
  // Backtrace and generate CIGAR
  banded_backtrace_matrix(banded_matrix,banded_pattern,text,text_length);
}

void banded_compute_unaligned(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded_unaligned(
      banded_matrix,banded_pattern,
      text,text_length);
  // Backtrace and generate CIGAR
  banded_backtrace_matrix_unaligned(banded_matrix,banded_pattern,text,text_length);
}

void banded_compute_blocking(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded_blocking(
      banded_matrix,banded_pattern,
      text,text_length);
  // Backtrace and generate CIGAR
  banded_backtrace_matrix_blocking(banded_matrix,banded_pattern,text,text_length);
}

void banded_compute_cutoff(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int text_length,
    const uint64_t cutoff_score) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded_cutoff(
      banded_matrix,banded_pattern,
      text,text_length,cutoff_score);
  // Backtrace and generate CIGAR
  banded_backtrace_matrix_cutoff(banded_matrix,banded_pattern,text,text_length,cutoff_score);
}

void banded_compute_cutoff_score(
    banded_matrix_t* const banded_matrix,
    banded_pattern_t* const banded_pattern,
    char* const text,
    const int64_t text_length,
    const int64_t text_finish_pos,
    const int64_t cutoff_score) {
  // Fill Matrix (Pv,Mv)
  bpm_compute_matrix_banded_cutoff_score(
      banded_matrix,banded_pattern,
      text,text_length,text_finish_pos,cutoff_score);
}

