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

#ifndef BPM_WINDOWED_H_
#define BPM_WINDOWED_H_

#include "utils/commons.h"
#include "alignment/cigar.h"
#include "system/mm_allocator.h"


typedef struct {
  /* BMP Pattern */
  char* pattern;                // Raw pattern
  uint64_t* PEQ;                // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_length;      // Length
  uint64_t pattern_num_words64; // ceil(Length / |w|)
  uint64_t pattern_mod;         // Length % |w|
  /* BPM Auxiliary data */
  uint64_t* P;
  uint64_t* M;
  uint64_t* level_mask;
  int64_t* score;
  int64_t* init_score;
  uint64_t* pattern_left;
} windowed_pattern_t;

/*
 * BPM matrix
 */
typedef struct {
  // Bit-encoded Matrix
  uint64_t* Pv;
  uint64_t* Mv;
  uint64_t min_score;
  uint64_t min_score_column;
  int64_t pos_v;
  int64_t pos_h;
  // CIGAR
  cigar_t* cigar;
} windowed_matrix_t;

/*
 * Setup
 */
void windowed_pattern_compile(
    windowed_pattern_t* const windowed_pattern,
    char* const pattern,
    const int pattern_length,
    mm_allocator_t* const mm_allocator);
void windowed_pattern_free(
    windowed_pattern_t* const windowed_pattern,
    mm_allocator_t* const mm_allocator);

void windowed_matrix_allocate(
    windowed_matrix_t* const windowed_matrix,
    const uint64_t pattern_length,
    const uint64_t text_length,
    mm_allocator_t* const mm_allocator,
    const int window_size);
void windowed_matrix_free(
    windowed_matrix_t* const windowed_matrix,
    mm_allocator_t* const mm_allocator);

/*
 * Edit distance computation using BPM
 */
void windowed_compute(
    windowed_matrix_t* const windowed_matrix,
    windowed_pattern_t* const windowed_pattern,
    char* const text,
    const int text_length,
    const int max_distance,
    const int window_size, 
    const int overlap_size);

#endif /* BPM_WINDOWED_H_ */

