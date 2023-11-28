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

void bpm_compute_matrix_hirschberg(
    const char* text,
    const char* text_r,
    const int64_t text_length,
    const char* pattern,
    const char* pattern_r,
    const int64_t pattern_length,
    const int64_t cutoff_score,
    cigar_t *cigar_out,
    mm_allocator_t *const mm_allocator)
{

    const int64_t k_end = ABS(((int64_t)text_length) - (int64_t)(pattern_length)) + 1;
    const int64_t cutoff_score_real = MAX(MAX(k_end, cutoff_score), 65);
    const int64_t sequence_length_diff = pattern_length - text_length;
    const int64_t relative_cutoff_score = DIV_CEIL((cutoff_score_real - ABS(sequence_length_diff)), 2);
    int64_t prolog_column_blocks;
    int64_t effective_bandwidth_blocks;
    if (sequence_length_diff >= 0)
    {
        prolog_column_blocks = DIV_CEIL(relative_cutoff_score, BPM_W64_LENGTH);
        effective_bandwidth_blocks = DIV_CEIL(relative_cutoff_score + sequence_length_diff, BPM_W64_LENGTH) + 1 + prolog_column_blocks;
    }
    else
    {
        prolog_column_blocks = DIV_CEIL(relative_cutoff_score - sequence_length_diff, BPM_W64_LENGTH);
        effective_bandwidth_blocks = DIV_CEIL(relative_cutoff_score, BPM_W64_LENGTH) + 1 + prolog_column_blocks;
    }

    uint64_t alignment_footprint = effective_bandwidth_blocks * text_length * BPM_W64_SIZE * 2;

    if (alignment_footprint > BUFFER_SIZE_16M)
    { // divide the alignment in 2

        const int64_t text_len = (text_length + 1) / 2;
        const int64_t text_len_r = text_length - text_len;

        const int64_t pattern_len = pattern_length;

        banded_pattern_t banded_pattern;
        banded_pattern_t banded_pattern_r;
        banded_pattern_compile(
            &banded_pattern, pattern,
            pattern_length, mm_allocator);
        banded_pattern_compile(
            &banded_pattern_r, pattern_r,
            pattern_length, mm_allocator);

        banded_matrix_t banded_matrix, banded_matrix_r;

        // Compute left side (for getting the central column)
        banded_matrix_allocate(
            &banded_matrix, pattern_length,
            text_length, cutoff_score, true, mm_allocator);

        banded_compute(
            &banded_matrix, &banded_pattern, text,
            text_length, text_len, true);

        // Compute right side (for getting the central column)
        banded_matrix_allocate(
            &banded_matrix_r, pattern_length,
            text_length, cutoff_score, true, mm_allocator);

        banded_compute(
            &banded_matrix_r, &banded_pattern_r, text_r,
            text_length, text_len_r, true);

        // vertival position of the first blocks computed on each aligments
        int64_t first_block_band_pos_v = text_len < prolog_column_blocks * BPM_W64_LENGTH ? 0 : (text_len / BPM_W64_LENGTH) - (prolog_column_blocks);
        int64_t first_block_band_pos_v_r = text_len_r < prolog_column_blocks * BPM_W64_LENGTH ? 0 : (text_len_r / BPM_W64_LENGTH) - (prolog_column_blocks);

        // Higher and lower cell's position computen in each aligments
        int64_t bottom_cell;
        int64_t higher_cell, higher_cell_r;
        int64_t starting_pos;
        const int64_t bottom_pos = banded_matrix.lower_block * 64 + 63 + first_block_band_pos_v * 64;
        const int64_t bottom_pos_r = (pattern_len - 1) - (banded_matrix_r.higher_block * 64 + 63 + first_block_band_pos_v_r * 64);
        const int64_t higher_pos = banded_matrix.higher_block * 64 + 63 + first_block_band_pos_v * 64;
        const int64_t higher_pos_r = (pattern_len - 1) - (banded_matrix_r.lower_block * 64 + 63 + first_block_band_pos_v_r * 64);

        // select lower cell between the two aligmnets
        if (bottom_pos > bottom_pos_r)
        {
            bottom_cell = banded_matrix.lower_block * 64 + 63;
            starting_pos = bottom_pos;
        }
        else
        {
            bottom_cell = bottom_pos_r - first_block_band_pos_v * 64;
            starting_pos = bottom_pos_r;
        }

        // select higher cell between the two aligmnets
        if (higher_pos < higher_pos_r)
        {
            higher_cell = banded_matrix.higher_block * 64 + 63;
            higher_cell_r = (pattern_len - 1) - higher_pos - first_block_band_pos_v_r * 64;
        }
        else
        {
            higher_cell = higher_pos_r - first_block_band_pos_v * 64;
            higher_cell_r = banded_matrix_r.lower_block * 64 + 63;
        }
        const uint64_t number_of_cells = higher_cell - bottom_cell + 2;

        int32_t *cell_score = (int32_t *)mm_allocator_malloc(mm_allocator, number_of_cells * sizeof(int32_t));
        int32_t *cell_score_r = (int32_t *)mm_allocator_malloc(mm_allocator, number_of_cells * sizeof(int32_t));

        // compute scores of the left side
        cell_score[0] = 0;
        for (uint64_t i = 0; i < number_of_cells; i++)
        {
            const uint64_t block = (bottom_cell + i) / BPM_W64_LENGTH;
            const uint64_t cell = (bottom_cell + i) % BPM_W64_LENGTH;
            cell_score[i + 1] = cell_score[i] + ((banded_matrix.Pv[block] >> cell) & 0x1ULL) - ((banded_matrix.Mv[block] >> cell) & 0x1ULL);
        }
        // compute scores of the right side
        cell_score_r[0] = 0;
        for (uint64_t i = 0; i < number_of_cells; i++)
        {
            const uint64_t block = (higher_cell_r + i) / BPM_W64_LENGTH;
            const uint64_t cell = (higher_cell_r + i) % BPM_W64_LENGTH;
            cell_score_r[i + 1] = cell_score_r[i] + ((banded_matrix_r.Pv[block] >> cell) & 0x1ULL) - ((banded_matrix_r.Mv[block] >> cell) & 0x1ULL);
        }

        // search the middle joint cell
        uint64_t smaller_pos = 0;
        int64_t smaller_score = cell_score_r[number_of_cells - 1] + cell_score[0];
        for (uint64_t i = 1; i < number_of_cells; i++)
        {
            int64_t new_score = cell_score_r[number_of_cells - 1 - i] + cell_score[i];
            if (new_score < smaller_score)
            {
                smaller_pos = i;
                smaller_score = new_score;
            }
        }

        // Divide the text and the pattern for the recursive call
        int64_t pattern_length_left = starting_pos + smaller_pos;
        int64_t pattern_length_right = pattern_length - pattern_length_left;

        const char* pattern_r_left = pattern_r + pattern_length_right;
        const char* pattern_right = pattern + pattern_length_left;

        int64_t text_length_right = text_length - text_len;
        const char* text_right = text + text_len;
        const char* text_r_left = text_r + text_length_right;

        // Obtain the score of rach sub aligmnet
        int64_t block_reference = DIV_CEIL(pattern_length_left, BPM_W64_LENGTH) - (number_of_cells < smaller_pos + BPM_W64_LENGTH);
        int64_t score_pos_l = block_reference * BPM_W64_LENGTH - (bottom_cell + first_block_band_pos_v * 64);
        int64_t score_l = cell_score[smaller_pos] - cell_score[score_pos_l] + banded_matrix.scores[block_reference - 1];

        int64_t block_reference_r = DIV_CEIL(pattern_length_right, BPM_W64_LENGTH) - (smaller_pos < BPM_W64_LENGTH);
        int64_t score_pos_r = block_reference_r * BPM_W64_LENGTH - (higher_cell_r + first_block_band_pos_v_r * 64);
        int64_t score_r = cell_score_r[number_of_cells - 1 - smaller_pos] - cell_score_r[score_pos_r] + banded_matrix_r.scores[block_reference_r - 1];

        // Free
        banded_pattern_free(&banded_pattern, mm_allocator);
        banded_pattern_free(&banded_pattern_r, mm_allocator);
        banded_matrix_free(&banded_matrix, mm_allocator);
        banded_matrix_free(&banded_matrix_r, mm_allocator);
        mm_allocator_free(mm_allocator, cell_score);
        mm_allocator_free(mm_allocator, cell_score_r);

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
    }
    else
    { // solve the alignment

        banded_pattern_t banded_pattern;
        banded_matrix_t banded_matrix;

        // Allocate
        banded_pattern_compile(
            &banded_pattern, pattern,
            pattern_length, mm_allocator);
        banded_matrix_allocate(
            &banded_matrix, pattern_length,
            text_length, cutoff_score, false, mm_allocator);

        // Align
        banded_compute(
            &banded_matrix, &banded_pattern, text,
            text_length, pattern_length, false);
        // Merge cigar
        cigar_prepend_forward(banded_matrix.cigar, cigar_out);
        // free variables
        banded_pattern_free(&banded_pattern, mm_allocator);
        banded_matrix_free(&banded_matrix, mm_allocator);
    }
}
