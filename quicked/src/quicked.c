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

#include "quicked.h"
#include "bpm_banded.h"
#include "bpm_windowed.h"
#include "bpm_hirschberg.h"
#include "utils/include/commons.h"
#include <stddef.h>

void extract_results(
    quicked_aligner_t *aligner,
    cigar_t *const cigar)
{
    if (aligner->params.only_score)
    {
        // Precomputed score
        aligner->score = cigar->score;
    }
    else
    {
        // CIGAR
        if (cigar->begin_offset < cigar->end_offset)
        {
            aligner->cigar = malloc((2 * (cigar->end_offset - cigar->begin_offset) + 10) * sizeof(char));
            cigar_sprint(aligner->cigar, cigar, true);
        }

        // Score from CIGAR
        aligner->score = cigar_score_edit(cigar);
    }
}

quicked_status_t run_banded(
    quicked_aligner_t *aligner,
    char *const pattern, const int pattern_len,
    char *const text, const int text_len)
{
    // FIXME: What if cutoff_score becomes 0?
    const int cutoff_score = (MAX(text_len, pattern_len) * (aligner->params.bandwidth)) / 100;

    // Allocate
    mm_allocator_t *const mm_allocator = aligner->mm_allocator;

    banded_pattern_t banded_pattern;
    banded_pattern_compile(&banded_pattern, pattern, pattern_len, mm_allocator);

    banded_matrix_t banded_matrix;
    banded_matrix_allocate(&banded_matrix, pattern_len, text_len, cutoff_score, aligner->params.only_score, mm_allocator);

    // Align
    // timer_start(&timer);
    banded_compute(&banded_matrix, &banded_pattern, text, text_len, text_len, aligner->params.only_score);
    // timer_stop(&timer);

    // Retrieve results
    extract_results(aligner, banded_matrix.cigar);

    // Free
    banded_pattern_free(&banded_pattern, mm_allocator);

    banded_matrix_free(&banded_matrix, mm_allocator);

    return QUICKED_WIP;
}

quicked_status_t run_windowed(
    quicked_aligner_t *aligner,
    char *const pattern, const int pattern_len,
    char *const text, const int text_len)
{
    const int window_size = aligner->params.window_size;
    const int overlap_size = aligner->params.overlap_size;

    // Allocate
    mm_allocator_t *const mm_allocator = aligner->mm_allocator;

    windowed_pattern_t windowed_pattern;
    windowed_pattern_compile(&windowed_pattern, pattern, pattern_len, mm_allocator);

    windowed_matrix_t windowed_matrix;
    windowed_matrix_allocate(&windowed_matrix, pattern_len, text_len, mm_allocator, window_size);

    // Align
    // timer_start(&timer);
    windowed_compute(&windowed_matrix, &windowed_pattern, text, window_size, overlap_size,
                     aligner->params.only_score, aligner->params.force_scalar);
    // timer_stop(&timer);

    // Retrieve results
    extract_results(aligner, windowed_matrix.cigar);

    // Free
    windowed_pattern_free(&windowed_pattern, mm_allocator);
    windowed_matrix_free(&windowed_matrix, mm_allocator);

    return QUICKED_WIP;
}

quicked_status_t run_hirschberg(
    quicked_aligner_t *aligner,
    char *const pattern, const int pattern_len,
    char *const text, const int text_len)
{
    // FIXME: What if cutoff_score becomes 0?
    const int cutoff_score = (MAX(text_len, pattern_len) * (aligner->params.bandwidth)) / 100;

    // Allocate
    mm_allocator_t *const mm_allocator = aligner->mm_allocator;

    char *text_r = (char *)mm_allocator_malloc(mm_allocator, text_len);
    char *pattern_r = (char *)mm_allocator_malloc(mm_allocator, pattern_len);

    reverse_string(text, text_r, text_len);
    reverse_string(pattern, pattern_r, pattern_len);

    cigar_t cigar_out;
    cigar_out.operations = (char *)malloc((pattern_len + text_len) * sizeof(char));
    cigar_out.begin_offset = pattern_len + text_len;
    cigar_out.end_offset = pattern_len + text_len;

    // Align
    // timer_start(&timer);
    bpm_compute_matrix_hirschberg(text, text_r, text_len, pattern, pattern_r, pattern_len,
                                  cutoff_score, &cigar_out, mm_allocator);
    // timer_stop(&timer);

    // Retrieve results
    extract_results(aligner, &cigar_out);

    // Free
    free(cigar_out.operations);
    return QUICKED_WIP;
}

quicked_status_t run_quicked(
    quicked_aligner_t *aligner,
    char *const pattern, const int pattern_len,
    char *const text, const int text_len)
{
    // TODO: Comment phases of the algorithm

    mm_allocator_t *const mm_allocator = aligner->mm_allocator;

    char *text_r = (char *)mm_allocator_malloc(mm_allocator, text_len);
    char *pattern_r = (char *)mm_allocator_malloc(mm_allocator, pattern_len);

    reverse_string(text, text_r, text_len);
    reverse_string(pattern, pattern_r, pattern_len);

    windowed_pattern_t windowed_pattern;
    windowed_pattern_compile(&windowed_pattern, pattern, pattern_len, mm_allocator);

    windowed_matrix_t windowed_matrix;
    windowed_matrix_allocate(&windowed_matrix, pattern_len, text_len, mm_allocator, 2);

    // timer_start(&align_input->timer);
    // timer_start(&align_input->timer_window_sse);

    // Align
    windowed_compute(&windowed_matrix, &windowed_pattern, text, 2, 1,
                     aligner->params.only_score, aligner->params.force_scalar);

    // timer_stop(&align_input->timer_window_sse);

    int64_t score = windowed_matrix.cigar->score;

    // Free
    windowed_pattern_free(&windowed_pattern, mm_allocator);
    windowed_matrix_free(&windowed_matrix, mm_allocator);

    if (score > MAX(text_len, pattern_len) / 4)
    {
        // timer_start(&align_input->timer_window_6x2);

        windowed_pattern_compile(&windowed_pattern, pattern, pattern_len, mm_allocator);
        windowed_matrix_allocate(&windowed_matrix, pattern_len, text_len, mm_allocator, 6);

        windowed_compute(&windowed_matrix, &windowed_pattern, text, 6, 1,
                         aligner->params.only_score, aligner->params.force_scalar);

        // align_input->seq_with_6x2 = true; // TODO: Remove if unused

        score = windowed_matrix.cigar->score;
        windowed_pattern_free(&windowed_pattern, mm_allocator);
        windowed_matrix_free(&windowed_matrix, mm_allocator);

        windowed_pattern_compile(&windowed_pattern, pattern_r, pattern_len, mm_allocator);
        windowed_matrix_allocate(&windowed_matrix, pattern_len, text_len, mm_allocator, 6);

        windowed_compute(&windowed_matrix, &windowed_pattern, text, 6, 1,
                         aligner->params.only_score, aligner->params.force_scalar);

        // align_input->seq_with_6x2_r = true; // TODO: Remove if unused

        score = MIN(score, windowed_matrix.cigar->score);

        windowed_pattern_free(&windowed_pattern, mm_allocator);
        windowed_matrix_free(&windowed_matrix, mm_allocator);

        // timer_stop(&align_input->timer_window_6x2);

        if (score > MAX(text_len, pattern_len) / 4)
        {
            // timer_start(&align_input->timer_banded_15);

            banded_pattern_t banded_pattern;
            banded_matrix_t banded_matrix_score;
            banded_pattern_compile(&banded_pattern, pattern, pattern_len, mm_allocator);

            score = MAX(text_len, pattern_len) * 3 / 20;

            banded_matrix_allocate(&banded_matrix_score, pattern_len, text_len, score, true, mm_allocator);

            banded_compute(&banded_matrix_score, &banded_pattern, text, text_len, text_len, true);

            // align_input->seqs_with_15 = true; // TODO: Remove if unused

            int64_t new_score = banded_matrix_score.cigar->score;

            banded_matrix_free(&banded_matrix_score, mm_allocator);

            // timer_stop(&align_input->timer_banded_15);

            while ((new_score > MAX(text_len, pattern_len) / 4 && score < new_score) || new_score < 0)
            {
                score *= 2;
                // timer_start(&align_input->timer_banded_30);

                banded_matrix_allocate(&banded_matrix_score, pattern_len, text_len, score, true, mm_allocator);

                banded_compute(&banded_matrix_score, &banded_pattern, text, text_len, text_len, true);

                // align_input->seqs_with_30 = true; // TODO: Remove if unused

                new_score = banded_matrix_score.cigar->score;

                banded_matrix_free(&banded_matrix_score, mm_allocator);

                // timer_stop(&align_input->timer_banded_30);
            }

            score = new_score;
        }
    }

    // timer_start(&align_input->timer_banded_hirschberg);

    cigar_t cigar_out;
    cigar_out.operations = (char *)malloc((pattern_len + text_len) * sizeof(char));
    cigar_out.begin_offset = pattern_len + text_len;
    cigar_out.end_offset = pattern_len + text_len;

    bpm_compute_matrix_hirschberg(text, text_r, text_len, pattern, pattern_r, pattern_len,
                                  score, &cigar_out, mm_allocator);

    // timer_stop(&align_input->timer_banded_hirschberg);
    // timer_stop(&align_input->timer);

    // benchmark_print_output(align_input, false, &cigar_out);
    // align_input->diff_scores = (float)(score - cigar_out.score) / (float)(MAX(align_input->text_length, align_input->pattern_length));

    extract_results(aligner, &cigar_out);

    free(cigar_out.operations);

    return QUICKED_WIP;
}

#pragma region API

quicked_params_t quicked_default_params()
{
    return (quicked_params_t){
        .algo = QUICKED,
        .only_score = false,
        .bandwidth = 15,
        .window_size = 2,
        .overlap_size = 1,
        .force_scalar = false,
    };
}

quicked_status_t quicked_new(
    quicked_aligner_t *aligner,
    quicked_params_t params)
{
    aligner->params = params;
    aligner->score = -1;
    aligner->cigar = NULL;
    aligner->mm_allocator = mm_allocator_new(BUFFER_SIZE_128M);

    return QUICKED_WIP;
}

quicked_status_t quicked_free(
    quicked_aligner_t *aligner)
{
    if (aligner->mm_allocator != NULL)
    {
        mm_allocator_delete(aligner->mm_allocator);
        aligner->mm_allocator = NULL;
    }

    if (aligner->cigar != NULL)
    {
        free(aligner->cigar);
        aligner->cigar = NULL;
    }

    return QUICKED_WIP;
}

quicked_status_t quicked_align(
    quicked_aligner_t *aligner,
    const char *pattern, const int pattern_len,
    const char *text, const int text_len)
{
    quicked_status_t status = QUICKED_ERROR;

    switch (aligner->params.algo)
    {
    case QUICKED:
        status = run_quicked(aligner, pattern, pattern_len, text, text_len);
        break;
    case WINDOWED:
        status = run_windowed(aligner, pattern, pattern_len, text, text_len);
        break;
    case BANDED:
        status = run_banded(aligner, pattern, pattern_len, text, text_len);
        break;
    case HIRSCHBERG:
        status = run_hirschberg(aligner, pattern, pattern_len, text, text_len);
        break;
    default:
        return QUICKED_UNKNOWN_ALGO;
    }

    return status;
}

#pragma endregion API