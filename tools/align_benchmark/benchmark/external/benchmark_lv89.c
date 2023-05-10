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

#include "benchmark/benchmark_utils.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_check.h"

#include "external/lv89/lv89.c"
#include "external/lv89/lv89.h"

/*
 * Benchmark LV89
 */
void benchmark_lv89(
    align_input_t* const align_input) {
  // Allocate
  uint8_t* mem = (uint8_t*)malloc((align_input->pattern_length+align_input->text_length)*16);
  // Align
  timer_start(&align_input->timer);
  const int score = lv_ed(
      align_input->pattern_length,align_input->pattern,
      align_input->text_length,align_input->text,
      false,mem);
  timer_stop(&align_input->timer);
  // Free
  free(mem);
  // Output. NOTE: No CIGAR is produced, just score
  cigar_t cigar = {.begin_offset=0,.end_offset=0,.score=score};
  if (align_input->output_file) {
    benchmark_print_output(align_input,true,&cigar);
  }
}
