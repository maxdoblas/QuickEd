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

#include "quicked_utils/include/commons.h"
#include "quicked_utils/include/mm_allocator.h"
#include "quicked_utils/include/dna_text.h"
#include "bpm_windowed.h"
#include "bpm_commons.h"

#ifndef __SSE4_1__
#warning "SSE4.1 or higher is required for SIMD windows"
#else
#include <immintrin.h>
#endif

/*
 * Setup
 */

void windowed_pattern_compile(
    windowed_pattern_t *const windowed_pattern,
    const char* pattern,
    const uint64_t pattern_length,
    mm_allocator_t *const mm_allocator)
{
    // Calculate dimensions
    const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length, BPM_W64_LENGTH);
    const uint64_t PEQ_length = pattern_num_words64 * BPM_W64_LENGTH;
    const uint64_t pattern_mod = pattern_length % BPM_W64_LENGTH;
    // Init fields
    windowed_pattern->pattern = pattern;
    windowed_pattern->pattern_length = pattern_length;
    windowed_pattern->pattern_num_words64 = pattern_num_words64;
    windowed_pattern->pattern_mod = pattern_mod;
    // Allocate memory
    const uint64_t aux_vector_size = pattern_num_words64 * BPM_W64_SIZE;
    const uint64_t PEQ_size = BPM_ALPHABET_LENGTH * aux_vector_size;
    const uint64_t score_size = pattern_num_words64 * UINT64_SIZE;
    const uint64_t total_memory = PEQ_size + 3 * aux_vector_size + 2 * score_size + (pattern_num_words64 + 1) * UINT64_SIZE;
    void *memory = mm_allocator_malloc(mm_allocator, total_memory);
    windowed_pattern->PEQ = memory;
    memory += PEQ_size;
    windowed_pattern->P = memory;
    memory += aux_vector_size;
    windowed_pattern->M = memory;
    memory += aux_vector_size;
    windowed_pattern->level_mask = memory;
    memory += aux_vector_size;
    windowed_pattern->score = memory;
    memory += score_size;
    windowed_pattern->init_score = memory;
    memory += score_size;
    windowed_pattern->pattern_left = memory;
    // Init PEQ
    memset(windowed_pattern->PEQ, 0, PEQ_size);
    uint64_t i;
    for (i = 0; i < pattern_length; ++i)
    {
        const uint8_t enc_char = dna_encode(pattern[i]);
        const uint64_t block = i / BPM_W64_LENGTH;
        const uint64_t mask = 1ull << (i % BPM_W64_LENGTH);
        windowed_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block, enc_char)] |= mask;
    }
    for (; i < PEQ_length; ++i)
    { // Padding
        const uint64_t block = i / BPM_W64_LENGTH;
        const uint64_t mask = 1ull << (i % BPM_W64_LENGTH);
        uint64_t j;
        for (j = 0; j < BPM_ALPHABET_LENGTH; ++j)
        {
            windowed_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block, j)] |= mask;
        }
    }
    // Init auxiliary data
    uint64_t pattern_left = pattern_length;
    const uint64_t top = pattern_num_words64 - 1;
    memset(windowed_pattern->level_mask, 0, aux_vector_size);
    for (i = 0; i < top; ++i)
    {
        windowed_pattern->level_mask[i] = BPM_W64_MASK;
        windowed_pattern->init_score[i] = BPM_W64_LENGTH;
        windowed_pattern->pattern_left[i] = pattern_left;
        pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left - BPM_W64_LENGTH : 0;
    }
    for (; i <= pattern_num_words64; ++i)
    {
        windowed_pattern->pattern_left[i] = pattern_left;
        pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left - BPM_W64_LENGTH : 0;
    }
    if (pattern_mod > 0)
    {
        const uint64_t mask_shift = pattern_mod - 1;
        windowed_pattern->level_mask[top] = 1ull << (mask_shift);
        windowed_pattern->init_score[top] = pattern_mod;
    }
    else
    {
        windowed_pattern->level_mask[top] = BPM_W64_MASK;
        windowed_pattern->init_score[top] = BPM_W64_LENGTH;
    }
}

void windowed_pattern_free(
    windowed_pattern_t *const windowed_pattern,
    mm_allocator_t *const mm_allocator)
{
    mm_allocator_free(mm_allocator, windowed_pattern->PEQ);
}

void windowed_matrix_allocate(
    windowed_matrix_t *const windowed_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    mm_allocator_t *const mm_allocator,
    const int window_size)
{
    // Parameters
    // const uint64_t num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
    const uint64_t num_words64 = window_size;
    // Allocate auxiliary matrix
    // const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
    const uint64_t aux_matrix_size = num_words64 * UINT64_SIZE * (BPM_W64_LENGTH * window_size + 2); /* (+1 base-column) */
    uint64_t *const Pv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
    uint64_t *const Mv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
    windowed_matrix->Mv = Mv;
    windowed_matrix->Pv = Pv;
    windowed_matrix->pos_v = pattern_length - 1;
    windowed_matrix->pos_h = text_length - 1;
    windowed_matrix->high_error_window = 0;
    // CIGAR
    windowed_matrix->cigar = cigar_new(pattern_length + text_length,mm_allocator);
    windowed_matrix->cigar->end_offset = pattern_length + text_length;
    windowed_matrix->cigar->begin_offset = pattern_length + text_length - 1;
    windowed_matrix->cigar->score = 0;

    const uint64_t aux_PEQ_size = window_size * UINT64_SIZE * BPM_ALPHABET_LENGTH; /* (+1 base-column) */
    windowed_matrix->PEQ_window = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_PEQ_size);
}

void windowed_matrix_free(
    windowed_matrix_t *const windowed_matrix,
    mm_allocator_t *const mm_allocator)
{
    mm_allocator_free(mm_allocator, windowed_matrix->Mv);
    mm_allocator_free(mm_allocator, windowed_matrix->Pv);
    mm_allocator_free(mm_allocator, windowed_matrix->PEQ_window);
    // CIGAR
    cigar_free(windowed_matrix->cigar,mm_allocator);
}

/*
 * Edit distance computation using BPM
 */

void windowed_reset_differences(
    uint64_t *const P,
    uint64_t *const M,
    const uint64_t max_distance)
{
    // Reset score,P,M
    for (uint64_t i = 0; i < max_distance; ++i)
    {
        P[i] = BPM_W64_ONES;
        M[i] = 0;
    }
}

void windowed_reset_differences_zero(
    uint64_t *const P,
    uint64_t *const M,
    const uint64_t max_distance)
{
    // Reset score,P,M
    for (uint64_t i = 0; i < max_distance; ++i)
    {
        P[i] = 0; // BPM_W64_ONES;
        M[i] = 0;
    }
}

void windowed_compute_window(
    windowed_matrix_t *const windowed_matrix,
    windowed_pattern_t *const windowed_pattern,
    const char* text,
    const int window_size)
{
    // Pattern variables
    const uint64_t *PEQ = windowed_pattern->PEQ;
    // const uint64_t num_words64 = windowed_pattern->pattern_num_words64;
    const uint64_t num_words64 = window_size;
    // int64_t* const score = windowed_pattern->score;
    uint64_t *const Pv = windowed_matrix->Pv;
    uint64_t *const Mv = windowed_matrix->Mv;
    uint64_t *const PEQ_window = windowed_matrix->PEQ_window;
    // const uint64_t max_distance__1 = max_distance+1;
    // Advance in DP-bit_encoded matrix
    int64_t text_position;
    int64_t pos_v_fi = windowed_matrix->pos_v;
    int64_t pos_h_fi = windowed_matrix->pos_h;

    int64_t pos_v = (pos_v_fi - UINT64_LENGTH * (window_size) + 1 >= 0) ? pos_v_fi - UINT64_LENGTH * (window_size) + 1 : 0;
    int64_t pos_h = (pos_h_fi - UINT64_LENGTH * (window_size) + 1 >= 0) ? pos_h_fi - UINT64_LENGTH * (window_size) + 1 : 0;

    if (pos_h == 0) {
        windowed_reset_differences(Pv, Mv, window_size);
    } else {
        windowed_reset_differences_zero(Pv, Mv, window_size);
    }

    int64_t steps_v = (pos_v_fi - pos_v) / UINT64_LENGTH + 1;
    int64_t steps_h = pos_h_fi - pos_h;
    uint64_t shift = pos_v % UINT64_LENGTH;
    uint64_t shift_mask = shift ? 0xFFFFFFFFFFFFFFFFULL : 0ULL;
    int64_t pos_v_block = (pos_v / UINT64_LENGTH);

    for (int64_t i = 0; i < steps_v; ++i)
    {
        for (uint64_t enc_char = 0; enc_char < BPM_ALPHABET_LENGTH; enc_char++)
        {
            const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i + pos_v_block, enc_char)] >> shift | ((PEQ[BPM_PATTERN_PEQ_IDX(i + pos_v_block + 1, enc_char)] << (BPM_W64_LENGTH - shift)) & shift_mask);
            PEQ_window[BPM_PATTERN_PEQ_IDX(i, enc_char)] = Eq;
        }
    }

    // First cell
    uint64_t Ph_first;
    if (pos_v == 0) {
        Ph_first = 1;
    } else {
        Ph_first = 0;
    }

    for (text_position = 0; text_position <= steps_h; ++text_position)
    {
        // Fetch next character
        const uint8_t enc_char = dna_encode(text[text_position + pos_h]);
        // Advance all blocks
        int64_t i;
        uint64_t PHin = Ph_first, MHin = 0, PHout, MHout;
        for (i = 0; i < steps_v; ++i)
        {
            /* Calculate Step Data */
            const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position, num_words64, i);
            const uint64_t next_bdp_idx = bdp_idx + num_words64;
            uint64_t Pv_in = Pv[bdp_idx];
            uint64_t Mv_in = Mv[bdp_idx];
            // const uint64_t mask = level_mask[i+pos_v_block];
            const uint64_t Eq = PEQ_window[BPM_PATTERN_PEQ_IDX(i, enc_char)];
            /* Compute Block */
            BPM_ADVANCE_BLOCK_NO_MASK(Eq, Pv_in, Mv_in, PHin, MHin, PHout, MHout);

            /* Adjust score and swap propagate Hv */
            Pv[next_bdp_idx] = Pv_in;
            Mv[next_bdp_idx] = Mv_in;
            PHin = PHout;
            MHin = MHout;
        }
    }
}

#ifdef __SSE4_1__
void windowed_compute_window_sse(
    windowed_matrix_t *const windowed_matrix,
    windowed_pattern_t *const windowed_pattern,
    const char* text,
    const int window_size)
{
    // Pattern variables
    const uint64_t *PEQ = windowed_pattern->PEQ;
    const uint64_t num_words64 = window_size;
    uint64_t *const Pv = windowed_matrix->Pv;
    uint64_t *const Mv = windowed_matrix->Mv;
    uint64_t *const PEQ_window = windowed_matrix->PEQ_window;
    // Advance in DP-bit_encoded matrix
    int64_t text_position;
    int64_t pos_v_fi = windowed_matrix->pos_v;
    int64_t pos_h_fi = windowed_matrix->pos_h;

    int64_t pos_v = (pos_v_fi - UINT64_LENGTH * (window_size) + 1 >= 0) ? pos_v_fi - UINT64_LENGTH * (window_size) + 1 : 0;
    int64_t pos_h = (pos_h_fi - UINT64_LENGTH * (window_size) + 1 >= 0) ? pos_h_fi - UINT64_LENGTH * (window_size) + 1 : 0;

    if (pos_h == 0) {
        windowed_reset_differences(Pv, Mv, window_size);
    } else {
        windowed_reset_differences_zero(Pv, Mv, window_size);
    }

    int64_t steps_v = (pos_v_fi - pos_v) / UINT64_LENGTH + 1;
    int64_t steps_h = pos_h_fi - pos_h;
    uint64_t shift = pos_v % UINT64_LENGTH;
    uint64_t shift_mask = shift ? 0xFFFFFFFFFFFFFFFFULL : 0ULL;
    int64_t pos_v_block = (pos_v / UINT64_LENGTH);

    // Generate aligned PEQ vectors
    for (int64_t i = 0; i < steps_v; ++i)
    {
        for (uint64_t enc_char = 0; enc_char < BPM_ALPHABET_LENGTH; enc_char++)
        {
            const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i + pos_v_block, enc_char)] >> shift | ((PEQ[BPM_PATTERN_PEQ_IDX(i + pos_v_block + 1, enc_char)] << (BPM_W64_LENGTH - shift)) & shift_mask);
            PEQ_window[BPM_PATTERN_PEQ_IDX(i, enc_char)] = Eq;
        }
    }

    // First cell
    uint64_t Ph_first, Mh_first = 0;
    if (pos_v == 0) {
        Ph_first = 1;
    } else {
        Ph_first = 0;
    }

    {
        const uint8_t enc_char = dna_encode(text[0 + pos_h]);
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(0, num_words64, 0);
        const uint64_t next_bdp_idx = bdp_idx + num_words64;
        uint64_t Pv_in = Pv[bdp_idx];
        uint64_t Mv_in = Mv[bdp_idx];
        uint64_t PHout, MHout;
        const uint64_t Eq = PEQ_window[BPM_PATTERN_PEQ_IDX(0, enc_char)];
        BPM_ADVANCE_BLOCK_NO_MASK(Eq, Pv_in, Mv_in, Ph_first, Mh_first, PHout, MHout);
        Pv[next_bdp_idx] = Pv_in;
        Mv[next_bdp_idx] = Mv_in;
        Ph_first = PHout;
        Mh_first = MHout;
    }
    __m128i PHout, MHout, PHin, MHin;
    PHin = _mm_set_epi64x(1, Ph_first);
    MHin = _mm_set_epi64x(0, Mh_first);

    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(0, num_words64, 1);
    uint64_t next_bdp_idx = bdp_idx + num_words64;
    __m128i Mv_v = _mm_loadu_si128((__m128i *)&Mv[bdp_idx]);
    __m128i Pv_v = _mm_loadu_si128((__m128i *)&Pv[bdp_idx]);

    for (text_position = 1; text_position <= steps_h; text_position += 2)
    {
        // Fetch next character
        uint8_t enc_char = dna_encode(text[text_position + pos_h]);
        uint8_t enc_char_2 = dna_encode(text[text_position + pos_h - 1]);
        uint8_t enc_char_3 = dna_encode(text[text_position + pos_h + 1]);

        /* Calculate Step Data */
        uint64_t Eq = PEQ_window[BPM_PATTERN_PEQ_IDX(0, enc_char)];
        uint64_t Eq_2 = PEQ_window[BPM_PATTERN_PEQ_IDX(1, enc_char_2)];

        /* Compute Block */
        __m128i Eq_v = _mm_set_epi64x(Eq, Eq_2);

        __m128i Xv = _mm_or_si128(Eq_v, Mv_v);
        __m128i _Eq = _mm_or_si128(Eq_v, MHin);
        __m128i sum_128 = _mm_add_epi64(_mm_and_si128(_Eq, Pv_v), Pv_v);
        __m128i Xh = _mm_or_si128(_mm_xor_si128(sum_128, Pv_v), _Eq);
        /* Calculate Hout */
        __m128i Ph_v = _mm_or_si128(Mv_v, ~(_mm_or_si128(Xh, Pv_v)));
        __m128i Mh_v = _mm_and_si128(Pv_v, Xh);
        /* Account Hout that propagates for the next block */
        PHout = _mm_srli_epi64(Ph_v, 63);
        MHout = _mm_srli_epi64(Mh_v, 63);
        /* Hout become the Hin of the next cell */
        Ph_v = _mm_slli_epi64(Ph_v, 1);
        Mh_v = _mm_slli_epi64(Mh_v, 1);
        /* Account Hin coming from the previous block */
        Ph_v = _mm_or_si128(Ph_v, PHin);
        Mh_v = _mm_or_si128(Mh_v, MHin);
        /* Finally, generate the Vout */
        Pv_v = _mm_or_si128(Mh_v, ~(_mm_or_si128(Xv, Ph_v)));
        Mv_v = _mm_and_si128(Ph_v, Xv);

        _mm_storeu_si128((__m128i *)&Pv[next_bdp_idx], Pv_v);
        _mm_storeu_si128((__m128i *)&Mv[next_bdp_idx], Mv_v);
        next_bdp_idx += num_words64;
        PHin = _mm_set_epi64x(1ULL, _mm_extract_epi64(PHout, 1));
        MHin = _mm_srli_si128(MHout, 8);

        // Second word (Unroll 2)
        Eq = PEQ_window[BPM_PATTERN_PEQ_IDX(0, enc_char_3)];
        Eq_2 = PEQ_window[BPM_PATTERN_PEQ_IDX(1, enc_char)];

        Eq_v = _mm_set_epi64x(Eq, Eq_2);

        Xv = _mm_or_si128(Eq_v, Mv_v);
        _Eq = _mm_or_si128(Eq_v, MHin);
        sum_128 = _mm_add_epi64(_mm_and_si128(_Eq, Pv_v), Pv_v);
        Xh = _mm_or_si128(_mm_xor_si128(sum_128, Pv_v), _Eq);
        /* Calculate Hout */
        Ph_v = _mm_or_si128(Mv_v, ~(_mm_or_si128(Xh, Pv_v)));
        Mh_v = _mm_and_si128(Pv_v, Xh);
        /* Account Hout that propagates for the next block */
        PHout = _mm_srli_epi64(Ph_v, 63);
        MHout = _mm_srli_epi64(Mh_v, 63);
        /* Hout become the Hin of the next cell */
        Ph_v = _mm_slli_epi64(Ph_v, 1);
        Mh_v = _mm_slli_epi64(Mh_v, 1);
        /* Account Hin coming from the previous block */
        Ph_v = _mm_or_si128(Ph_v, PHin);
        Mh_v = _mm_or_si128(Mh_v, MHin);
        /* Finally, generate the Vout */
        Pv_v = _mm_or_si128(Mh_v, ~(_mm_or_si128(Xv, Ph_v)));
        Mv_v = _mm_and_si128(Ph_v, Xv);

        _mm_storeu_si128((__m128i *)&Pv[next_bdp_idx], Pv_v);
        _mm_storeu_si128((__m128i *)&Mv[next_bdp_idx], Mv_v);
        PHin = _mm_set_epi64x(0ULL, _mm_extract_epi64(PHout, 1));
        MHin = _mm_srli_si128(MHout, 8);
        next_bdp_idx += num_words64;
    }
    // Last cell
    Ph_first = _mm_extract_epi64(PHout, 1);
    Mh_first = _mm_extract_epi64(MHout, 1);
    {
        const uint8_t enc_char = dna_encode(text[steps_h + pos_h]);
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(steps_h, num_words64, 1);
        const uint64_t next_bdp_idx = bdp_idx + num_words64;
        uint64_t Pv_in = Pv[bdp_idx];
        uint64_t Mv_in = Mv[bdp_idx];
        uint64_t PHout, MHout;
        const uint64_t Eq = PEQ_window[BPM_PATTERN_PEQ_IDX(1, enc_char)];
        BPM_ADVANCE_BLOCK_NO_MASK(Eq, Pv_in, Mv_in, Ph_first, Mh_first, PHout, MHout);
        Pv[next_bdp_idx] = Pv_in;
        Mv[next_bdp_idx] = Mv_in;
        Ph_first = PHout;
        Mh_first = MHout;
    }
}
#endif

void windowed_backtrace(
    windowed_matrix_t *const windowed_matrix,
    const windowed_pattern_t *const windowed_pattern,
    const char* text,
    const int window_size,
    const int overlap_size)
{
    // Parameters
    const char* pattern = windowed_pattern->pattern;
    const uint64_t *const Pv = windowed_matrix->Pv;
    const uint64_t *const Mv = windowed_matrix->Mv;
    char *const operations = windowed_matrix->cigar->operations;
    int op_sentinel = windowed_matrix->cigar->begin_offset;
    // Retrieve the alignment. Store the match
    const uint64_t num_words64 = window_size;
    int64_t h = windowed_matrix->pos_h;
    int64_t v = windowed_matrix->pos_v;
    int64_t h_min = windowed_matrix->pos_h - UINT64_LENGTH * (window_size) + 1 > 0 ? (windowed_matrix->pos_h - (window_size)*UINT64_LENGTH + 1) : 0;
    int64_t h_overlap = windowed_matrix->pos_h - UINT64_LENGTH * (window_size - overlap_size) + 1 > 0 ? (windowed_matrix->pos_h - (window_size - overlap_size) * UINT64_LENGTH + 1) : 0;
    int64_t v_min = windowed_matrix->pos_v - UINT64_LENGTH * (window_size) + 1 > 0 ? (windowed_matrix->pos_v - (window_size)*UINT64_LENGTH + 1) : 0;
    int64_t v_overlap = windowed_matrix->pos_v - UINT64_LENGTH * (window_size - overlap_size) + 1 > 0 ? (windowed_matrix->pos_v - (window_size - overlap_size) * UINT64_LENGTH + 1) : 0;

    while (v >= v_overlap && h >= h_overlap)
    {
        const uint8_t block = (v - v_min) / UINT64_LENGTH;
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX((h - h_min + 1), num_words64, block);
        const uint64_t mask = 1L << (v - v_min % UINT64_LENGTH);

        if (text[h] == pattern[v])
        {
            operations[op_sentinel--] = 'M';
            --h;
            --v;
        } else if (Pv[bdp_idx] & mask)
        {
            operations[op_sentinel--] = 'D';
            --v;
        }
        else if (Mv[(bdp_idx - num_words64)] & mask)
        {
            operations[op_sentinel--] = 'I';
            --h;
        }
        else
        {
            operations[op_sentinel--] = 'X';
            --h;
            --v;
        }
    }
    windowed_matrix->pos_h = h;
    windowed_matrix->pos_v = v;

    windowed_matrix->cigar->begin_offset = op_sentinel;
}

void windowed_backtrace_score_only(
    windowed_matrix_t *const windowed_matrix,
    const windowed_pattern_t *const windowed_pattern,
    const char* text,
    const int hew_threshold,
    const int window_size,
    const int overlap_size)
{
    // Parameters
    const char* pattern = windowed_pattern->pattern;
    const uint64_t *const Pv = windowed_matrix->Pv;
    const uint64_t *const Mv = windowed_matrix->Mv;
    // Retrieve the alignment. Store the match
    const uint64_t num_words64 = window_size;
    int64_t h = windowed_matrix->pos_h;
    int64_t v = windowed_matrix->pos_v;
    int64_t h_min = windowed_matrix->pos_h - UINT64_LENGTH * (window_size) + 1 > 0 ? (windowed_matrix->pos_h - (window_size)*UINT64_LENGTH + 1) : 0;
    int64_t h_overlap = windowed_matrix->pos_h - UINT64_LENGTH * (window_size - overlap_size) + 1 > 0 ? (windowed_matrix->pos_h - (window_size - overlap_size) * UINT64_LENGTH + 1) : 0;
    int64_t v_min = windowed_matrix->pos_v - UINT64_LENGTH * (window_size) + 1 > 0 ? (windowed_matrix->pos_v - (window_size)*UINT64_LENGTH + 1) : 0;
    int64_t v_overlap = windowed_matrix->pos_v - UINT64_LENGTH * (window_size - overlap_size) + 1 > 0 ? (windowed_matrix->pos_v - (window_size - overlap_size) * UINT64_LENGTH + 1) : 0;
    int64_t score = 0;

    while (v >= v_overlap && h >= h_overlap)
    {
        const uint8_t block = (v - v_min) / UINT64_LENGTH;
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX((h - h_min + 1), num_words64, block);
        const uint64_t mask = 1L << (v - v_min % UINT64_LENGTH);

        if (Pv[bdp_idx] & mask)
        {
            score++;
            --v;
        }
        else if (Mv[(bdp_idx - num_words64)] & mask)
        {
            score++;
            --h;
        }
        else if ((text[h] == pattern[v]))
        {
            --h;
            --v;
        }
        else
        {
            score++;
            --h;
            --v;
        }
    }

    if (score > ((window_size - overlap_size) * UINT64_LENGTH * hew_threshold / 100))
        windowed_matrix->high_error_window++;

    windowed_matrix->pos_h = h;
    windowed_matrix->pos_v = v;
    windowed_matrix->cigar->score += score;
}

void windowed_compute(
    windowed_matrix_t *const windowed_matrix,
    windowed_pattern_t *const windowed_pattern,
    const char* text,
    const int hew_threshold,
    const int window_size,
    const int overlap_size,
    const bool score_only,
    const bool force_scalar)
{
    while (windowed_matrix->pos_v >= 0 && windowed_matrix->pos_h >= 0)
    {
        // Fill window (Pv,Mv)
        #ifdef __SSE4_1__
        if (!force_scalar && window_size == 2) // Vectorized version only works for window_size == 2
        {
            windowed_compute_window_sse(windowed_matrix, windowed_pattern, text, window_size);
        }
        else
        #endif
        {
            windowed_compute_window(windowed_matrix, windowed_pattern, text, window_size);
        }

        // Compute window backtrace
        if (score_only)
        {
            windowed_backtrace_score_only(windowed_matrix, windowed_pattern, text, hew_threshold, window_size, overlap_size);
        }
        else
        {
            windowed_backtrace(windowed_matrix, windowed_pattern, text, window_size, overlap_size);
        }
    }

    if (score_only)
    {
        int64_t h = windowed_matrix->pos_h;
        int64_t v = windowed_matrix->pos_v;
        if (h >= 0)
            windowed_matrix->cigar->score += h + 1;
        if (v >= 0)
            windowed_matrix->cigar->score += v + 1;
    }
    else
    {
        int64_t h = windowed_matrix->pos_h;
        int64_t v = windowed_matrix->pos_v;
        char *const operations = windowed_matrix->cigar->operations;
        int op_sentinel = windowed_matrix->cigar->begin_offset;
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
        windowed_matrix->pos_h = h;
        windowed_matrix->pos_v = v;
        windowed_matrix->cigar->begin_offset = op_sentinel + 1;
    }
}
