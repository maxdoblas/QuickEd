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

#ifndef QUICKED_H
#define QUICKED_H

#include <quicked_utils/include/mm_allocator.h>
#include <quicked_utils/include/profiler_timer.h>
#include <stdbool.h>

#define QUICKED_WINDOW_STAGES 2 // Number of window sizes to go through before doing banded
#define QUICKED_FAST_WINDOW_SIZE 2
#define QUICKED_FAST_WINDOW_OVERLAP 1

typedef enum {
    QUICKED,
    WINDOWED,
    BANDED,
    HIRSCHBERG,
} quicked_algo_t;

typedef struct quicked_params_t {
    quicked_algo_t algo;
    unsigned int bandwidth;
    unsigned int window_size;
    unsigned int overlap_size;
    unsigned int hew_threshold[QUICKED_WINDOW_STAGES];
    unsigned int hew_percentage[QUICKED_WINDOW_STAGES];
    bool only_score;
    bool force_scalar;
    bool external_timer;
    mm_allocator_t *external_allocator;
} quicked_params_t;

typedef struct quicked_aligner_t {
    quicked_params_t* params;
    mm_allocator_t *mm_allocator;
    char* cigar;
    int score;
    // Profiling
    profiler_timer_t *timer;
    profiler_timer_t *timer_windowed_s;
    profiler_timer_t *timer_windowed_l;
    profiler_timer_t *timer_banded;
    profiler_timer_t *timer_align;
} quicked_aligner_t;

typedef enum quicked_status_t {
    QUICKED_OK                   = 0,
    QUICKED_ERROR                = -1,  // Default error code
    QUICKED_FAIL_NON_CONVERGENCE = -2,  // The hirschberg has no solution for the actual cutoff score
    QUICKED_UNKNOWN_ALGO         = -3,  // Provided algorithm is not supported
    QUICKED_EMPTY_SEQUENCE       = -4,  // Empty sequence

    // Development codes
    QUICKED_UNIMPLEMENTED        = -10, // Function declared but not implemented
    QUICKED_WIP                  = 1,   // Function implementation in progress. Considered not an error
} quicked_status_t;

bool quicked_check_error(quicked_status_t status);
const char* quicked_status_msg(quicked_status_t status);

quicked_params_t quicked_default_params(void);
quicked_status_t quicked_new(
    quicked_aligner_t *aligner,
    quicked_params_t *params
);
quicked_status_t quicked_free(
    quicked_aligner_t *aligner
);
quicked_status_t quicked_align(
    quicked_aligner_t *aligner,
    const char* pattern, const int pattern_len,
    const char* text, const int text_len
);

#endif // QUICKED_H
