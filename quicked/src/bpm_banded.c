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

#include "quicked_utils/include/commons.h"
#include "quicked_utils/include/mm_allocator.h"
#include "quicked_utils/include/dna_text.h"
#include "bpm_banded.h"
#include "bpm_commons.h"
#include <sys/mman.h>
#include <immintrin.h>

#define BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv, Mv, PHin, MHin, PHout, MHout) 	       \
    __m256i Xv    = _mm256_or_si256(Eq, Mv);      /*Eq | Mv*/                          \
    __m256i _Eq   = _mm256_or_si256(Eq, MHin);    /*Eq | MHin*/      	       	       \
    __m256i Xh    = _mm256_and_si256(_Eq, Pv);    /*(((_Eq & Pv) + Pv) ^ Pv) | _Eq*/   \
    	    Xh    = _mm256_add_epi64(Xh, Pv);					                       \
    	    Xh    = _mm256_xor_si256(Xh, Pv); 					                       \
    	    Xh    = _mm256_or_si256(Xh, _Eq);					                       \
    __m256i Ph    = _mm256_or_si256(Xh, Pv);      /*Mv | ~(Xh | Pv)*/                  \
    	    Ph    = _mm256_or_si256(Mv, ~Ph);					                       \
    __m256i Mh 	  = _mm256_and_si256(Pv, Xh);     /*Pv & Xh*/			               \
            PHout = _mm256_and_si256(Ph, mask);   /*(Ph & mask) != 0*/	        	   \
    	    PHout = _mm256_cmpeq_epi64(PHout, _mm256_setzero_si256());		           \
            PHout = _mm256_andnot_si256(PHout, _mm256_set1_epi64x(1));				   \
            MHout = _mm256_and_si256(Mh, mask);   /*(Mh & mask) != 0 */	               \
    	    MHout = _mm256_cmpeq_epi64(MHout, _mm256_setzero_si256());				   \
            MHout = _mm256_andnot_si256(MHout, _mm256_set1_epi64x(1));				   \
      	    Ph 	  = _mm256_slli_epi64(Ph, 1);     /*Ph <<= 1*/		            	   \
      	    Mh 	  = _mm256_slli_epi64(Mh, 1);     /*Mh <<= 1*/		            	   \
      	    Ph    = _mm256_or_si256(Ph, PHin);    /*Ph |= PHin*/	        		   \
      	    Mh    = _mm256_or_si256(Mh, MHin);    /*Mh |= MHin*/		        	   \
      	    Pv    = _mm256_or_si256(Xv, Ph);      /*Mh | ~(Xv | Ph)*/      		       \
      	    Pv    = _mm256_or_si256(Mh, ~Pv);                                          \
            Mv    = _mm256_and_si256(Ph, Xv);                                          \


#define BPM_ADVANCE_BLOCK_SI256_2(Eq2, mask2, Pv2, Mv2, PHin2, MHin2, PHout2, MHout2) 	       \
    __m256i Xv2    = _mm256_or_si256(Eq2, Mv2);      /*Eq | Mv*/                          \
    __m256i _Eq2   = _mm256_or_si256(Eq2, MHin2);    /*Eq | MHin*/      	       	       \
    __m256i Xh2    = _mm256_and_si256(_Eq2, Pv2);    /*(((_Eq & Pv) + Pv) ^ Pv) | _Eq*/   \
    	    Xh2    = _mm256_add_epi64(Xh2, Pv2);					                       \
    	    Xh2    = _mm256_xor_si256(Xh2, Pv2); 					                       \
    	    Xh2    = _mm256_or_si256(Xh2, _Eq2);					                       \
    __m256i Ph2    = _mm256_or_si256(Xh2, Pv2);      /*Mv | ~(Xh | Pv)*/                  \
    	    Ph2    = _mm256_or_si256(Mv2, ~Ph2);					                       \
    __m256i Mh2 	  = _mm256_and_si256(Pv2, Xh2);     /*Pv & Xh*/			               \
            PHout2 = _mm256_and_si256(Ph2, mask2);   /*(Ph & mask) != 0*/	        	   \
    	    PHout2 = _mm256_cmpeq_epi64(PHout2, _mm256_setzero_si256());		           \
            PHout2 = _mm256_andnot_si256(PHout2, _mm256_set1_epi64x(1));				   \
            MHout2 = _mm256_and_si256(Mh2, mask2);   /*(Mh & mask) != 0 */	               \
    	    MHout2 = _mm256_cmpeq_epi64(MHout2, _mm256_setzero_si256());				   \
            MHout2 = _mm256_andnot_si256(MHout2, _mm256_set1_epi64x(1));				   \
      	    Ph2 	  = _mm256_slli_epi64(Ph2, 1);     /*Ph <<= 1*/		            	   \
      	    Mh2 	  = _mm256_slli_epi64(Mh2, 1);     /*Mh <<= 1*/		            	   \
      	    Ph2    = _mm256_or_si256(Ph2, PHin2);    /*Ph |= PHin*/	        		   \
      	    Mh2    = _mm256_or_si256(Mh2, MHin2);    /*Mh |= MHin*/		        	   \
      	    Pv2    = _mm256_or_si256(Xv2, Ph2);      /*Mh | ~(Xv | Ph)*/      		       \
      	    Pv2    = _mm256_or_si256(Mh2, ~Pv2);                                          \
            Mv2    = _mm256_and_si256(Ph2, Xv2);                                          \

void banded_pattern_compile(
    banded_pattern_t *const banded_pattern,
    const char* pattern,
    const uint64_t pattern_length,
    mm_allocator_t *const mm_allocator)
{
    // Calculate dimensions
    const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length, BPM_W64_LENGTH);
    const uint64_t PEQ_length = pattern_num_words64 * BPM_W64_LENGTH;
    const uint64_t pattern_mod = pattern_length % BPM_W64_LENGTH;
    // Init fields
    banded_pattern->pattern = pattern;
    banded_pattern->pattern_length = pattern_length;
    banded_pattern->pattern_num_words64 = pattern_num_words64;
    banded_pattern->pattern_mod = pattern_mod;
    // Allocate memory
    const uint64_t aux_vector_size = pattern_num_words64 * BPM_W64_SIZE;
    const uint64_t PEQ_size = BPM_ALPHABET_LENGTH * aux_vector_size;
    const uint64_t total_memory = PEQ_size + 3 * aux_vector_size + (pattern_num_words64)*UINT64_SIZE;
    void *memory = mm_allocator_malloc(mm_allocator, total_memory);
    banded_pattern->PEQ = memory;
    memory += PEQ_size;
    banded_pattern->P = memory;
    memory += aux_vector_size;
    banded_pattern->M = memory;
    memory += aux_vector_size;
    banded_pattern->level_mask = memory;
    // Init PEQ
    memset(banded_pattern->PEQ, 0, PEQ_size);
    uint64_t i;
    for (i = 0; i < pattern_length; ++i)
    {
        const uint8_t enc_char = dna_encode(pattern[i]);
        const uint64_t block = i / BPM_W64_LENGTH;
        const uint64_t mask = 1ull << (i % BPM_W64_LENGTH);
        banded_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block, enc_char)] |= mask;
    }
    for (; i < PEQ_length; ++i)
    { // Padding
        const uint64_t block = i / BPM_W64_LENGTH;
        const uint64_t mask = 1ull << (i % BPM_W64_LENGTH);
        uint64_t j;
        for (j = 0; j < BPM_ALPHABET_LENGTH; ++j)
        {
            banded_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block, j)] |= mask;
        }
    }
    // Init auxiliary data
    const uint64_t top = pattern_num_words64 - 1;
    memset(banded_pattern->level_mask, 0, aux_vector_size);
    for (i = 0; i < top; ++i)
    {
        banded_pattern->level_mask[i] = BPM_W64_MASK;
    }
    if (pattern_mod > 0)
    {
        const uint64_t mask_shift = pattern_mod - 1;
        banded_pattern->level_mask[top] = 1ull << (mask_shift);
    }
    else
    {
        banded_pattern->level_mask[top] = BPM_W64_MASK;
    }
}

void banded_pattern_free(
    banded_pattern_t *const banded_pattern,
    mm_allocator_t *const mm_allocator)
{
    mm_allocator_free(mm_allocator, banded_pattern->PEQ);
}

void banded_matrix_allocate(
    banded_matrix_t *const banded_matrix,
    const int64_t pattern_length,
    const int64_t text_length,
    const int64_t cutoff_score,
    bool only_score,
    mm_allocator_t *const mm_allocator)
{

    const int64_t k_end = ABS(((int64_t)text_length) - (int64_t)(pattern_length)) + 1;
    banded_matrix->cutoff_score = MAX(MAX(k_end, cutoff_score), 65);
    banded_matrix->sequence_length_diff = pattern_length - text_length;
    banded_matrix->relative_cutoff_score = DIV_CEIL((banded_matrix->cutoff_score - ABS(banded_matrix->sequence_length_diff)), 2);
    if (banded_matrix->sequence_length_diff >= 0)
    {
        banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH);
        banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score + banded_matrix->sequence_length_diff, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
    }
    else
    {
        banded_matrix->prolog_column_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score - banded_matrix->sequence_length_diff, BPM_W64_LENGTH);
        banded_matrix->effective_bandwidth_blocks = DIV_CEIL(banded_matrix->relative_cutoff_score, BPM_W64_LENGTH) + 1 + banded_matrix->prolog_column_blocks;
    }
    banded_matrix->effective_bandwidth = banded_matrix->cutoff_score;

    const uint64_t num_words64 = banded_matrix->effective_bandwidth_blocks;
    // Allocate auxiliary matrix
    const uint64_t num_cols = (only_score ? 1 : (text_length + 1)); // Only 1 column if only score, or text_length + 1 columns
    const uint64_t aux_matrix_size = num_words64 * UINT64_SIZE * num_cols;
    uint64_t * Pv;
    uint64_t * Mv;
    if (aux_matrix_size > BUFFER_SIZE_256K){
        Pv = (uint64_t *)mm_allocator_allocate(mm_allocator, aux_matrix_size,false,BUFFER_SIZE_2M);
        Mv = (uint64_t *)mm_allocator_allocate(mm_allocator, aux_matrix_size,false,BUFFER_SIZE_2M);
        #ifdef __linux__
            if (madvise(Pv, aux_matrix_size, MADV_HUGEPAGE) == -1) {
                perror("madvise");
                exit(EXIT_FAILURE);
            }
            if (madvise(Mv, aux_matrix_size, MADV_HUGEPAGE) == -1) {
                perror("madvise");
                exit(EXIT_FAILURE);
            }
        #endif
    } else {
        Pv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
        Mv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
    }
    int64_t *const scores = (int64_t *)mm_allocator_malloc(mm_allocator, (DIV_CEIL(pattern_length, BPM_W64_LENGTH) + num_words64 / 2) * UINT64_SIZE);
    banded_matrix->Mv = Mv;
    banded_matrix->Pv = Pv;
    banded_matrix->scores = scores;
    // CIGAR
    banded_matrix->cigar = cigar_new(pattern_length + text_length,mm_allocator);
    banded_matrix->cigar->end_offset = pattern_length + text_length;
}

void banded_matrix_free(
    banded_matrix_t *const banded_matrix,
    mm_allocator_t *const mm_allocator)
{
    mm_allocator_free(mm_allocator, banded_matrix->Mv);
    mm_allocator_free(mm_allocator, banded_matrix->Pv);
    mm_allocator_free(mm_allocator, banded_matrix->scores);
    // CIGAR
    cigar_free(banded_matrix->cigar,mm_allocator);
}

void bpm_reset_search(
    const uint64_t num_words,
    uint64_t *const P,
    uint64_t *const M,
    int64_t *const scores)
{
    // Reset P,M
    uint64_t i;
    P[0] = BPM_W64_ONES;
    M[0] = 0;
    scores[0] = BPM_W64_LENGTH;
    for (i = 1; i < num_words; ++i)
    {
        P[i] = BPM_W64_ONES;
        M[i] = 0;
        scores[i]  = scores[i-1] + BPM_W64_LENGTH;
    }
}

void bpm_compute_matrix_banded_cutoff(
    banded_matrix_t *const banded_matrix,
    banded_pattern_t *const banded_pattern,
    const char* text,
    const int64_t text_length)
{
    // Pattern variables
    const uint64_t *PEQ = banded_pattern->PEQ;
    const int64_t effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;

    const int64_t num_block_rows = DIV_CEIL(banded_pattern->pattern_length, BPM_W64_LENGTH);

    const uint64_t *const level_mask = banded_pattern->level_mask;
    uint64_t *const Pv = banded_matrix->Pv;
    uint64_t *const Mv = banded_matrix->Mv;
    int64_t *const scores = banded_matrix->scores;
    const uint64_t num_words64 = effective_bandwidth_blocks;

    const int64_t sequence_length_diff = banded_matrix->sequence_length_diff;
    const int64_t prologue_columns = banded_matrix->prolog_column_blocks;
    const int64_t finish_v_pos_inside_band = prologue_columns * BPM_W64_LENGTH + sequence_length_diff;

    // Prepare last block of the next column
    int64_t pos_v = -prologue_columns;
    int64_t pos_h = 0;
    int64_t first_block_v = prologue_columns;
    int64_t last_block_v = effective_bandwidth_blocks - 1;

    bpm_reset_search(effective_bandwidth_blocks, Pv, Mv, scores);

    // Advance in DP-bit_encoded matrix
    int64_t text_position;
    // Main loop
    for (text_position = 0; text_position < text_length; ++text_position)
    {
        // Fetch next character
        const uint8_t enc_char = dna_encode(text[text_position]);
        // Advance all blocks
        int64_t i;
        uint64_t PHin = 1, MHin = 0, PHout, MHout;
        // Main Loop
        for (i = first_block_v; i <= last_block_v; ++i)
        {
 
            /* Calculate Step Data */
            const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position, num_words64, i);
            const uint64_t next_bdp_idx = bdp_idx + num_words64;
            uint64_t Pv_in = Pv[bdp_idx];
            uint64_t Mv_in = Mv[bdp_idx];
            const uint64_t mask = level_mask[i + pos_v];
            const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v), enc_char)];

            /* Compute Block */
            BPM_ADVANCE_BLOCK(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

            /* Swap propagate Hv */
            Pv[next_bdp_idx] = Pv_in;
            Mv[next_bdp_idx] = Mv_in;

            PHin = PHout;
            MHin = MHout;
            scores[i + pos_v] = scores[i + pos_v] + PHout - MHout;
        }

        // Update the score for the new column
        if ((text_position + 1) % 64 == 0)
        {
            // printf("-----------------------------------------------------\n");
            //  chech if the band of the lower side should be cutted
            int cut_band_lower = (first_block_v + 2 < last_block_v) && (finish_v_pos_inside_band > BPM_W64_LENGTH * (first_block_v + 1)) && (scores[first_block_v + pos_v + 1] + (finish_v_pos_inside_band - BPM_W64_LENGTH * (first_block_v + 1))) > banded_matrix->cutoff_score;

            if (cut_band_lower && (pos_h >= prologue_columns))
            {
                first_block_v++;
            }
            else if (!cut_band_lower && (pos_h < prologue_columns))
            {
                first_block_v--;
            }

            uint64_t next_bdp_idx = BPM_PATTERN_BDP_IDX(text_position + 1, num_words64, 0);
            // Shift results one block in the last column of a 64-column block
            for (int64_t j = first_block_v; j < last_block_v; j++)
            {
                Pv[next_bdp_idx + j] = Pv[next_bdp_idx + j + 1];
                Mv[next_bdp_idx + j] = Mv[next_bdp_idx + j + 1];
            }
            Pv[next_bdp_idx + last_block_v] = BPM_W64_ONES;
            Mv[next_bdp_idx + last_block_v] = 0;
            // Update the score for the new column
            int64_t pos = last_block_v + pos_v;
            scores[pos + 1] = scores[pos] + BPM_W64_LENGTH;

            // chech if the band of the higher side should be cutted
            int cut_band_higer = (first_block_v + 2 < last_block_v) && (BPM_W64_LENGTH * (last_block_v - 1) > finish_v_pos_inside_band) && (scores[last_block_v + pos_v - 1] + (BPM_W64_LENGTH * (last_block_v - 1) - finish_v_pos_inside_band)) > banded_matrix->cutoff_score;

            if (cut_band_higer || (pos_v + last_block_v >= num_block_rows - 1))
            {
                last_block_v--;
            }
            pos_v++;
            pos_h++;
        }
    }

    uint64_t final_score;
    if (banded_pattern->pattern_length % BPM_W64_LENGTH)
    {
        final_score = scores[(banded_pattern->pattern_length) / BPM_W64_LENGTH] - (BPM_W64_LENGTH - (banded_pattern->pattern_length % BPM_W64_LENGTH));
    }
    else
    {
        final_score = scores[(banded_pattern->pattern_length - 1) / BPM_W64_LENGTH];
    }
    banded_matrix->cigar->score = final_score;
    banded_matrix->higher_block = last_block_v;
    banded_matrix->lower_block = first_block_v;
}

void print(uint64_t* v, int64_t start, int64_t end, int k) {
    puts("PRINT PV");
    printf("ITER: %d\n", k);
    for(int64_t i = start; i <= end; i++) {
        printf("Pv[%ld] = %ld\n", i, v[i]);
    }
}

static inline __attribute__((always_inline)) void compute_advance_block (
    uint64_t* Pv, 
    uint64_t* Mv,
    const uint64_t *const PEQ,
    const uint64_t *const level_mask,
    int64_t* scores,
    uint64_t i,
    uint64_t pos_v, 
    uint8_t enc_char, 
    uint64_t* PHin, 
    uint64_t* MHin
) 
{
    uint64_t PHout, MHout; 
    uint64_t Pv_in = Pv[i];
    uint64_t Mv_in = Mv[i];
    uint64_t mask  = level_mask[i + pos_v];
    uint64_t Eq    = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v), enc_char)];
    uint64_t _PHin = *PHin; 
    uint64_t _MHin = *MHin;
    
    BPM_ADVANCE_BLOCK(Eq, mask, Pv_in, Mv_in, _PHin, _MHin, PHout, MHout);

    Pv[i] = Pv_in;
    Mv[i] = Mv_in;
    *PHin = PHout;
    *MHin = MHout;
    scores[i + pos_v] = scores[i + pos_v] + PHout - MHout;
}

void bpm_compute_matrix_banded_cutoff_score(
    banded_matrix_t *const banded_matrix,
    banded_pattern_t *const banded_pattern,
    const char* text,
    const int64_t text_length,
    const int64_t text_finish_pos)
{
    // Pattern variables
    const uint64_t *PEQ = banded_pattern->PEQ;
    // TODO: remove if necessary
    const int64_t k_end = ABS(((int64_t)text_length) - (int64_t)(banded_pattern->pattern_length)) + 1;
    const int64_t real_bandwidth = MAX(MAX(k_end, banded_matrix->cutoff_score), 65);
    const int64_t effective_bandwidth_blocks = DIV_CEIL(real_bandwidth, BPM_W64_LENGTH) + 1;
    const int64_t num_block_rows = DIV_CEIL(banded_pattern->pattern_length, BPM_W64_LENGTH);

    const uint64_t *const level_mask = banded_pattern->level_mask;
    uint64_t *const Pv = banded_matrix->Pv;
    uint64_t *const Mv = banded_matrix->Mv;
    int64_t *const scores = banded_matrix->scores;

    const int64_t sequence_length_diff = banded_matrix->sequence_length_diff;
    const int64_t prologue_columns = banded_matrix->prolog_column_blocks;
    const int64_t finish_v_pos_inside_band = prologue_columns * BPM_W64_LENGTH + sequence_length_diff;

    // Prepare last block of the next column
    int64_t pos_v = -prologue_columns;
    int64_t pos_h = 0;
    int64_t first_block_v = prologue_columns;
    int64_t last_block_v = effective_bandwidth_blocks - 1;

    bpm_reset_search(effective_bandwidth_blocks, Pv, Mv, scores);

    // Advance in DP-bit_encoded matrix
    int64_t text_position, k;
    // uint64_t count = 0;
    //  Main loop

    int64_t text_block = (text_finish_pos / 64);

    for (k = 0; k < text_block; k++)
    {
       if (last_block_v - first_block_v >= 8)
        {
            for (text_position = k * 64; text_position < (k+1) * 64; text_position+=8)
            {
                // Fetch next character
                const uint8_t enc_char1 = dna_encode(text[text_position]);
                const uint8_t enc_char2 = dna_encode(text[text_position+1]);
                const uint8_t enc_char3 = dna_encode(text[text_position+2]);
                const uint8_t enc_char4 = dna_encode(text[text_position+3]);
                const uint8_t enc_char5 = dna_encode(text[text_position+4]);
                const uint8_t enc_char6 = dna_encode(text[text_position+5]);
                const uint8_t enc_char7 = dna_encode(text[text_position+6]);
                const uint8_t enc_char8 = dna_encode(text[text_position+7]);
                

                int64_t i = first_block_v;
                
                uint64_t PHin_0 = 1, MHin_0 = 0;
                uint64_t PHin_1 = 1, MHin_1 = 0;
                uint64_t PHin_2 = 1, MHin_2 = 0; 
                uint64_t PHin_3 = 1, MHin_3 = 0;
                uint64_t PHin_4 = 1, MHin_4 = 0;
                uint64_t PHin_5 = 1, MHin_5 = 0;
                uint64_t PHin_6 = 1, MHin_6 = 0; 
                uint64_t PHin_7 = 1, MHin_7 = 0;

                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+2, pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char3, &PHin_2, &MHin_2);
                
                __m256i Pv_in = _mm256_set_epi64x(Pv[first_block_v+2], Pv[first_block_v+1], Pv[first_block_v], 0); 
                __m256i Mv_in = _mm256_set_epi64x(Mv[first_block_v+2], Mv[first_block_v+1], Mv[first_block_v], 0); 
                __m256i PHin  = _mm256_set_epi64x(PHin_0, PHin_1, PHin_2, PHin_3);  
                __m256i MHin  = _mm256_set_epi64x(MHin_0, MHin_1, MHin_2, MHin_3);
                __m256i MHout = MHin; 
                __m256i PHout = PHout;
                
                #pragma GCC unroll(4)
                for (i = first_block_v+3; i < first_block_v+7; i++) 
                {
                    Pv_in = _mm256_permute4x64_epi64(Pv_in, 0x39);
                    Mv_in = _mm256_permute4x64_epi64(Mv_in, 0x39);

                    Pv_in = _mm256_insert_epi64(Pv_in,   Pv[i], 3);
                    Mv_in = _mm256_insert_epi64(Mv_in,   Mv[i], 3);

                    uint64_t Eq_3 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v),     enc_char1)];
                    uint64_t Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 1), enc_char2)];
                    uint64_t Eq_1 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 2), enc_char3)];
                    uint64_t Eq_0 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 3), enc_char4)];
                    
                    __m256i Eq    = _mm256_set_epi64x (Eq_3, Eq_2, Eq_1, Eq_0); 
                    __m256i score = _mm256_lddqu_si256((__m256i*)&scores[i+pos_v-3]);
                    __m256i mask  = _mm256_lddqu_si256((__m256i const*)&level_mask[i+pos_v-3]);
                    
                    BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

                    Pv[i-3] = _mm256_extract_epi64(Pv_in, 0);
                    Mv[i-3] = _mm256_extract_epi64(Mv_in, 0);

                    MHin = MHout; 
                    PHin = PHout;
                    score = _mm256_add_epi64(score, PHout);
                    score = _mm256_sub_epi64(score, MHout);
                    _mm256_storeu_si256((__m256i*)&scores[i+pos_v-3], score);
               }
                _mm256_storeu_si256((__m256i*)&Pv[first_block_v+3], Pv_in);
                _mm256_storeu_si256((__m256i*)&Mv[first_block_v+3], Mv_in);

                PHin_0 = _mm256_extract_epi64(PHout, 3);
                MHin_0 = _mm256_extract_epi64(MHout, 3);
                PHin_1 = _mm256_extract_epi64(PHout, 2);
                MHin_1 = _mm256_extract_epi64(MHout, 2);
                PHin_2 = _mm256_extract_epi64(PHout, 1);
                MHin_2 = _mm256_extract_epi64(MHout, 1);
                PHin_3 = _mm256_extract_epi64(PHout, 0);
                MHin_3 = _mm256_extract_epi64(MHout, 0);

                i = first_block_v;
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char5, &PHin_4, &MHin_4);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char5, &PHin_4, &MHin_4);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+2, pos_v, enc_char5, &PHin_4, &MHin_4);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char6, &PHin_5, &MHin_5);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char6, &PHin_5, &MHin_5);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char7, &PHin_6, &MHin_6);

                PHin = _mm256_set_epi64x(PHin_4, PHin_5, PHin_6, PHin_7);
                MHin = _mm256_set_epi64x(MHin_4, MHin_5, MHin_6, MHin_7);
                Pv_in = _mm256_set_epi64x(Pv[first_block_v+2], Pv[first_block_v+1], Pv[first_block_v], 0); 
                Mv_in = _mm256_set_epi64x(Mv[first_block_v+2], Mv[first_block_v+1], Mv[first_block_v], 0); 
               
                __m256i PHin2 = PHout;
                __m256i MHin2 = MHout;
                __m256i PHout2 = PHout, MHout2 = MHout;
                __m256i Pv_in2 = _mm256_set_epi64x(Pv[first_block_v+6], Pv[first_block_v+5], Pv[first_block_v+4], 0); 
                __m256i Mv_in2 = _mm256_set_epi64x(Mv[first_block_v+6], Mv[first_block_v+5], Mv[first_block_v+4], 0); 

                // Main Loop
                for (i = first_block_v+7; i <= last_block_v; ++i)
                {                    
                    Pv_in2 = _mm256_permute4x64_epi64(Pv_in2, 0x39);
                    Mv_in2 = _mm256_permute4x64_epi64(Mv_in2, 0x39);

                    Pv_in2 = _mm256_insert_epi64(Pv_in2,   Pv[i], 3);
                    Mv_in2 = _mm256_insert_epi64(Mv_in2,   Mv[i], 3);

                    uint64_t Eq_3 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v),     enc_char1)];
                    uint64_t Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 1), enc_char2)];
                    uint64_t Eq_1 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 2), enc_char3)];
                    uint64_t Eq_0 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 3), enc_char4)];
                    __m256i Eq    = _mm256_set_epi64x (Eq_3, Eq_2, Eq_1, Eq_0); 
                    __m256i score = _mm256_lddqu_si256((__m256i*)&scores[i+pos_v-3]);
                    __m256i mask  = _mm256_lddqu_si256((__m256i const*)&level_mask[i+pos_v-3]);
                    
                    BPM_ADVANCE_BLOCK_SI256_2(Eq, mask, Pv_in2, Mv_in2, PHin2, MHin2, PHout2, MHout2);

                    Pv[i-3] = _mm256_extract_epi64(Pv_in2, 0);
                    Mv[i-3] = _mm256_extract_epi64(Mv_in2, 0);

                    MHin2 = MHout2; 
                    PHin2 = PHout2;
                    score = _mm256_add_epi64(score, PHout2);
                    score = _mm256_sub_epi64(score, MHout2);
                    _mm256_storeu_si256((__m256i*)&scores[i+pos_v-3], score);

                    Pv_in = _mm256_permute4x64_epi64(Pv_in, 0x39);
                    Mv_in = _mm256_permute4x64_epi64(Mv_in, 0x39);

                    Pv_in = _mm256_insert_epi64(Pv_in,   Pv[i-4], 3);
                    Mv_in = _mm256_insert_epi64(Mv_in,   Mv[i-4], 3);

                    Eq_3 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 4), enc_char5)];
                    Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 5), enc_char6)];
                    Eq_1 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 6), enc_char7)];
                    Eq_0 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 7), enc_char8)];
                    Eq    = _mm256_set_epi64x (Eq_3, Eq_2, Eq_1, Eq_0); 
                    score = _mm256_lddqu_si256((__m256i*)&scores[i+pos_v-7]);
                    mask  = _mm256_lddqu_si256((__m256i const*)&level_mask[i+pos_v-7]);
                    
                    BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

                    Pv[i-7] = _mm256_extract_epi64(Pv_in, 0);
                    Mv[i-7] = _mm256_extract_epi64(Mv_in, 0);

                    MHin = MHout; 
                    PHin = PHout;
                    score = _mm256_add_epi64(score, PHout);
                    score = _mm256_sub_epi64(score, MHout);
                    _mm256_storeu_si256((__m256i*)&scores[i+pos_v-7], score);
                }
                _mm256_storeu_si256((__m256i*)&Pv[last_block_v-3], Pv_in2);
                _mm256_storeu_si256((__m256i*)&Mv[last_block_v-3], Mv_in2);
                _mm256_storeu_si256((__m256i*)&Pv[last_block_v-7], Pv_in);
                _mm256_storeu_si256((__m256i*)&Mv[last_block_v-7], Mv_in);

                PHin_1 = _mm256_extract_epi64(PHout2, 2);
                MHin_1 = _mm256_extract_epi64(MHout2, 2);
                PHin_2 = _mm256_extract_epi64(PHout2, 1);
                MHin_2 = _mm256_extract_epi64(MHout2, 1);
                PHin_3 = _mm256_extract_epi64(PHout2, 0);
                MHin_3 = _mm256_extract_epi64(MHout2, 0);


                i = last_block_v;
                
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char3, &PHin_2, &MHin_2);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char3, &PHin_2, &MHin_2);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-2, pos_v, enc_char4, &PHin_3, &MHin_3);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char4, &PHin_3, &MHin_3);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char4, &PHin_3, &MHin_3);
                
                Pv_in = _mm256_set_epi64x(Pv[i-4], Pv[i-5], Pv[i-6], Pv[i-7]); 
                Mv_in = _mm256_set_epi64x(Mv[i-4], Mv[i-5], Mv[i-6], Mv[i-7]); 
                MHout = MHin; 
                PHout = PHin;

                #pragma GCC unroll(4)
                for (i = last_block_v-3; i <= last_block_v; i++) 
                {
                    Pv_in = _mm256_permute4x64_epi64(Pv_in, 0x39);
                    Mv_in = _mm256_permute4x64_epi64(Mv_in, 0x39);

                    Pv_in = _mm256_insert_epi64(Pv_in,   Pv[i], 3);
                    Mv_in = _mm256_insert_epi64(Mv_in,   Mv[i], 3);

                    uint64_t Eq_3 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v),     enc_char5)];
                    uint64_t Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 1), enc_char6)];
                    uint64_t Eq_1 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 2), enc_char7)];
                    uint64_t Eq_0 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 3), enc_char8)];
                    __m256i Eq    = _mm256_set_epi64x (Eq_3, Eq_2, Eq_1, Eq_0); 
                    __m256i score = _mm256_lddqu_si256((__m256i*)&scores[i+pos_v-3]);
                    __m256i mask  = _mm256_lddqu_si256((__m256i const*)&level_mask[i+pos_v-3]);
                    
                    BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

                    Pv[i-3] = _mm256_extract_epi64(Pv_in, 0);
                    Mv[i-3] = _mm256_extract_epi64(Mv_in, 0);

                    MHin = MHout; 
                    PHin = PHout;
                    score = _mm256_add_epi64(score, PHout);
                    score = _mm256_sub_epi64(score, MHout);
                    _mm256_storeu_si256((__m256i*)&scores[i+pos_v-3], score);
               }
                _mm256_storeu_si256((__m256i*)&Pv[last_block_v-3], Pv_in);
                _mm256_storeu_si256((__m256i*)&Mv[last_block_v-3], Mv_in);

                PHin_5 = _mm256_extract_epi64(PHout, 2);
                MHin_5 = _mm256_extract_epi64(MHout, 2);
                PHin_6 = _mm256_extract_epi64(PHout, 1);
                MHin_6 = _mm256_extract_epi64(MHout, 1);
                PHin_7 = _mm256_extract_epi64(PHout, 0);
                MHin_7 = _mm256_extract_epi64(MHout, 0);

                i = last_block_v;

                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char6, &PHin_5, &MHin_5);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char7, &PHin_6, &MHin_6);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char7, &PHin_6, &MHin_6);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-2, pos_v, enc_char8, &PHin_7, &MHin_7);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char8, &PHin_7, &MHin_7);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char8, &PHin_7, &MHin_7);
                }
        }
        else if (last_block_v - first_block_v >= 4)
        {
            for (text_position = k * 64; text_position < (k+1) * 64; text_position+=4)
            {
                // Fetch next character
                const uint8_t enc_char1 = dna_encode(text[text_position]);
                const uint8_t enc_char2 = dna_encode(text[text_position+1]);
                const uint8_t enc_char3 = dna_encode(text[text_position+2]);
                const uint8_t enc_char4 = dna_encode(text[text_position+3]);

                int64_t i = first_block_v;
                
                uint64_t PHin_0 = 1ul, MHin_0 = 0ul, PHin_1 = 1ul, MHin_1 = 0ul, PHin_2 = 1ul, MHin_2 = 0ul, PHin_3 = 1ul, MHin_3 = 0ul;

                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+2, pos_v, enc_char1, &PHin_0, &MHin_0);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i+1, pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char3, &PHin_2, &MHin_2);

                __m256i Pv_in = _mm256_set_epi64x(Pv[first_block_v+2], Pv[first_block_v+1], Pv[first_block_v], 0); 
                __m256i Mv_in = _mm256_set_epi64x(Mv[first_block_v+2], Mv[first_block_v+1], Mv[first_block_v], 0); 
                __m256i PHin  = _mm256_set_epi64x(PHin_0, PHin_1, PHin_2, PHin_3);  
                __m256i MHin  = _mm256_set_epi64x(MHin_0, MHin_1, MHin_2, MHin_3);
                __m256i MHout = MHin; 
                __m256i PHout = PHout;

                for (i = first_block_v+3; i <= last_block_v; ++i)
                {
                    Pv_in = _mm256_permute4x64_epi64(Pv_in, 0x39);
                    Mv_in = _mm256_permute4x64_epi64(Mv_in, 0x39);

                    Pv_in = _mm256_insert_epi64(Pv_in,   Pv[i], 3);
                    Mv_in = _mm256_insert_epi64(Mv_in,   Mv[i], 3);

                    uint64_t Eq_3 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v),     enc_char1)];
                    uint64_t Eq_2 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 1), enc_char2)];
                    uint64_t Eq_1 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 2), enc_char3)];
                    uint64_t Eq_0 = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v - 3), enc_char4)];
                    
                    __m256i Eq    = _mm256_set_epi64x (Eq_3, Eq_2, Eq_1, Eq_0); 
                    __m256i score = _mm256_lddqu_si256((__m256i*)&scores[i+pos_v-3]);
                    __m256i mask  = _mm256_lddqu_si256((__m256i const*)&level_mask[i+pos_v-3]);
                    
                    BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

                    Pv[i-3] = _mm256_extract_epi64(Pv_in, 0);
                    Mv[i-3] = _mm256_extract_epi64(Mv_in, 0);

                    MHin = MHout; 
                    PHin = PHout;
                    score = _mm256_add_epi64(score, PHout);
                    score = _mm256_sub_epi64(score, MHout);
                    _mm256_storeu_si256((__m256i*)&scores[i+pos_v-3], score);
                }
                i = last_block_v;
                _mm256_storeu_si256((__m256i*)&Pv[i-3], Pv_in);
                _mm256_storeu_si256((__m256i*)&Mv[i-3], Mv_in);

                PHin_1 = _mm256_extract_epi64(PHout, 2);
                MHin_1 = _mm256_extract_epi64(MHout, 2);
                PHin_2 = _mm256_extract_epi64(PHout, 1);
                MHin_2 = _mm256_extract_epi64(MHout, 1);
                PHin_3 = _mm256_extract_epi64(PHout, 0);
                MHin_3 = _mm256_extract_epi64(MHout, 0);

                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char2, &PHin_1, &MHin_1);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char3, &PHin_2, &MHin_2);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char3, &PHin_2, &MHin_2);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-2, pos_v, enc_char4, &PHin_3, &MHin_3);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char4, &PHin_3, &MHin_3);
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v, enc_char4, &PHin_3, &MHin_3);
            }          
        }
        else 
        {
            for (text_position = k * 64; text_position < (k+1) * 64; text_position+=2)
            {
                const uint8_t enc_char  = dna_encode(text[text_position]);
                const uint8_t enc_char2 = dna_encode(text[text_position+1]);
                int64_t i = first_block_v;  
                uint64_t PHin_0 = 1, MHin_0 = 0, PHin_1 = 1, MHin_1 = 0;
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i, pos_v, enc_char, &PHin_0, &MHin_0);
                for (i = first_block_v+1; i <= last_block_v; ++i)
                {
                    compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i-1, pos_v, enc_char2, &PHin_1, &MHin_1);
                    compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i,   pos_v,  enc_char, &PHin_0, &MHin_0);
                }
                i = last_block_v;
                compute_advance_block(Pv, Mv, PEQ, level_mask, scores, i, pos_v, enc_char2, &PHin_1, &MHin_1);
            }
        }
        

        // chech if the band of the lower side should be cutted
        int cut_band_lower = (first_block_v + 2 < last_block_v) && (finish_v_pos_inside_band > BPM_W64_LENGTH * (first_block_v + 1)) && (scores[first_block_v + pos_v + 1] + (finish_v_pos_inside_band - BPM_W64_LENGTH * (first_block_v + 1))) > banded_matrix->cutoff_score;

        // if we are in the prolog columns, we have to decrease the first_block_v (apart of cutting the band if necessary)
        if (cut_band_lower && (pos_h >= prologue_columns))
        {
            first_block_v++;
        }
        else if (!cut_band_lower && (pos_h < prologue_columns))
        {
            first_block_v--;
        }

        // Shift results one block in the last column of a 64-column block
        for (int64_t j = first_block_v; j < last_block_v; j++)
        {
            Pv[j] = Pv[j + 1];
            Mv[j] = Mv[j + 1];
        }
        Pv[last_block_v] = BPM_W64_ONES;
        Mv[last_block_v] = 0;
        // Update the score for the new column
        int64_t pos = last_block_v + pos_v;
        scores[pos + 1] = scores[pos] + BPM_W64_LENGTH;

        // chech if the band of the higher side should be cutted
        int cut_band_higer = (first_block_v + 2 < last_block_v) && (BPM_W64_LENGTH * (last_block_v - 1) > finish_v_pos_inside_band) && (scores[last_block_v + pos_v - 1] + (BPM_W64_LENGTH * (last_block_v - 1) - finish_v_pos_inside_band)) > banded_matrix->cutoff_score;

        if (cut_band_higer || (pos_v + last_block_v >= num_block_rows))
        {
            last_block_v--;
        }
        pos_v++;
        pos_h++;
    }

    for (; text_position < text_finish_pos; ++text_position)
    {
        // Fetch next character
        const uint8_t enc_char = dna_encode(text[text_position]);
        // Advance all blocks
        int64_t i;
        uint64_t PHin = 1, MHin = 0, PHout, MHout;
        // Main Loop
        for (i = first_block_v; i <= last_block_v; ++i)
        {
            /* Calculate Step Data */
            uint64_t Pv_in = Pv[i];
            uint64_t Mv_in = Mv[i];
            const uint64_t mask = level_mask[i + pos_v];
            const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX((i + pos_v), enc_char)];

            /* Compute Block */
            BPM_ADVANCE_BLOCK(Eq, mask, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

            /* Swap propagate Hv */
            Pv[i] = Pv_in;
            Mv[i] = Mv_in;
            PHin = PHout;
            MHin = MHout;
            scores[i + pos_v] = scores[i + pos_v] + PHout - MHout;
        }
    }
    uint64_t final_score;
    if (banded_pattern->pattern_length % BPM_W64_LENGTH)
    {
        final_score = scores[(banded_pattern->pattern_length) / BPM_W64_LENGTH] - (BPM_W64_LENGTH - (banded_pattern->pattern_length % BPM_W64_LENGTH));
    }
    else
    {
        final_score = scores[(banded_pattern->pattern_length - 1) / BPM_W64_LENGTH];
    }
    banded_matrix->cigar->score = final_score;
    banded_matrix->higher_block = last_block_v;
    banded_matrix->lower_block = first_block_v;
}

void banded_backtrace_matrix_cutoff(
    banded_matrix_t *const banded_matrix,
    const banded_pattern_t *const banded_pattern,
    const char* text,
    const int64_t text_length)
{
    // Parameters
    const char* pattern = banded_pattern->pattern;
    const uint64_t pattern_length = banded_pattern->pattern_length;
    const uint64_t *const Pv = banded_matrix->Pv;
    const uint64_t *const Mv = banded_matrix->Mv;
    char *const operations = banded_matrix->cigar->operations;
    int op_sentinel = banded_matrix->cigar->end_offset - 1;
    const int effective_bandwidth_blocks = banded_matrix->effective_bandwidth_blocks;
    const int64_t prologue_columns = banded_matrix->prolog_column_blocks;

    // Retrieve the alignment. Store the match
    const uint64_t num_words64 = effective_bandwidth_blocks;
    int64_t h = text_length - 1;
    int64_t v = pattern_length - 1;

    while (v >= 0 && h >= 0)
    {
        const int64_t block_h = h / BPM_W64_LENGTH;
        const int64_t block_h_r = (h + 1) / BPM_W64_LENGTH;
        const int64_t effective_v = v - BPM_W64_LENGTH * (block_h - prologue_columns);
        const int64_t effective_v_r = v - BPM_W64_LENGTH * (block_h_r - prologue_columns);
        const int64_t block_v = effective_v / BPM_W64_LENGTH;
        const int64_t block_v_r = effective_v_r / BPM_W64_LENGTH;
        const int64_t bdp_idx = BPM_PATTERN_BDP_IDX(h, num_words64, block_v);
        const int64_t bdp_idx_r = BPM_PATTERN_BDP_IDX(h + 1, num_words64, block_v_r);
        const uint64_t mask = 1UL << (effective_v % BPM_W64_LENGTH);
        const uint64_t mask_r = 1UL << (effective_v_r % BPM_W64_LENGTH);

        // CIGAR operation Test
        if (Pv[bdp_idx_r] & mask_r)
        {
            operations[op_sentinel--] = 'D';
            --v;
        }
        else if (Mv[(bdp_idx)] & mask)
        {
            operations[op_sentinel--] = 'I';
            --h;
        }
        else if ((text[h] == pattern[v]))
        {
            operations[op_sentinel--] = 'M';
            --h;
            --v;
        }
        else
        {
            operations[op_sentinel--] = 'X';
            --h;
            --v;
        }
    }
    while (h >= 0)
    {
        operations[op_sentinel--] = 'I';
        --h;
    }
    while (v >= 0)
    {
        operations[op_sentinel--] = 'D';
        --v;
    }
    banded_matrix->cigar->begin_offset = op_sentinel + 1;
}

void banded_compute(
    banded_matrix_t *const banded_matrix,
    banded_pattern_t *const banded_pattern,
    const char* text,
    const int64_t text_length,
    const int64_t text_finish_pos,
    const bool only_score)
{
    if (only_score)
    {
        // Fill Matrix (Pv,Mv)
        bpm_compute_matrix_banded_cutoff_score(banded_matrix, banded_pattern, text, text_length, text_finish_pos);
    }
    else
    {
        // Fill Matrix (Pv,Mv)
        bpm_compute_matrix_banded_cutoff(banded_matrix, banded_pattern, text, text_length);

        // Backtrace and generate CIGAR
        banded_backtrace_matrix_cutoff(banded_matrix, banded_pattern, text, text_length);
    }
}
