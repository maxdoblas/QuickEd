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

#include "utils/include/commons.h"
#include "utils/include/mm_allocator.h"
#include "utils/include/dna_text.h"
#include "bpm_banded.h"
#include "bpm_commons.h"

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

    const int64_t num_words64 = banded_matrix->effective_bandwidth_blocks;
    // Allocate auxiliary matrix
    const int64_t num_cols = (only_score ? 1 : (text_length + 1)); // Only 1 column if only score, or text_length + 1 columns
    const int64_t aux_matrix_size = num_words64 * UINT64_SIZE * num_cols;
    uint64_t *const Pv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
    uint64_t *const Mv = (uint64_t *)mm_allocator_malloc(mm_allocator, aux_matrix_size);
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
        scores[i] = scores[i - 1] + BPM_W64_LENGTH;
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
    int64_t text_position;
    // uint64_t count = 0;
    //  Main loop
    for (text_position = 0; text_position < text_finish_pos; ++text_position)
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

        // Update the score for the new column

        if ((text_position + 1) % 64 == 0)
        {

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