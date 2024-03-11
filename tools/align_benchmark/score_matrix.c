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
#include "score_matrix.h"

/*
 * Setup
 */
void score_matrix_allocate(
    score_matrix_t* const score_matrix,
    const int num_rows,
    const int num_columns,
    mm_allocator_t* const mm_allocator) {
  // Allocate DP matrix
  int h;
  score_matrix->num_rows = num_rows;
  score_matrix->num_columns = num_columns;
  score_matrix->columns = mm_allocator_malloc(mm_allocator,num_columns*sizeof(int*)); // Columns
  for (h=0;h<num_columns;++h) {
    score_matrix->columns[h] = mm_allocator_calloc(mm_allocator,num_rows,int,false); // Rows
  }
  // MM
  score_matrix->mm_allocator = mm_allocator;
}
void score_matrix_free(
    score_matrix_t* const score_matrix) {
  // Parameters
  mm_allocator_t* const mm_allocator = score_matrix->mm_allocator;
  // DP matrix
  const int num_columns = score_matrix->num_columns;
  int h;
  for (h=0;h<num_columns;++h) {
    mm_allocator_free(mm_allocator,score_matrix->columns[h]);
  }
  mm_allocator_free(mm_allocator,score_matrix->columns);
}
  
