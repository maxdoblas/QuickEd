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
    } else {
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

    if (aligner->params.only_score)
    {
        banded_matrix_allocate_cutoff_score(&banded_matrix, pattern_len, text_len, cutoff_score, mm_allocator);

        // Align
        // timer_start(&align_input->timer);
        banded_compute_cutoff_score(&banded_matrix, &banded_pattern, text, text_len, text_len);
        // timer_stop(&align_input->timer);
    }
    else
    {
        banded_matrix_allocate_cutoff(&banded_matrix, pattern_len, text_len, cutoff_score, mm_allocator);
        // Align
        // timer_start(&align_input->timer);
        banded_compute_cutoff(&banded_matrix, &banded_pattern, text, text_len);
        // timer_stop(&align_input->timer);
    }

    // Retrieve results
    extract_results(aligner, banded_matrix.cigar);

    // Free
    banded_pattern_free(&banded_pattern, mm_allocator);

    banded_matrix_free_cutoff(&banded_matrix, mm_allocator);

    return QUICKED_WIP;
}

quicked_params_t quicked_default_params()
{
    return (quicked_params_t){
        .algo = QUICKED,
        .only_score = false,
        .bandwidth = 0};
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
    char *const pattern, const int pattern_len,
    char *const text, const int text_len)
{
    quicked_status_t status = QUICKED_ERROR;

    switch (aligner->params.algo)
    {
    case QUICKED:
        break;
    case WINDOWED:
        break;
    case BANDED:
        status = run_banded(aligner, pattern, pattern_len, text, text_len);
        break;
    case ERA:
        break;
    case H_ERA:
        break;
    default:
        return QUICKED_UNKNOWN_ALGO;
    }

    return status;
}