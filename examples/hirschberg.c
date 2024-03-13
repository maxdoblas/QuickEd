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

#include <quicked.h>
#include <stdio.h>
#include <string.h>

int main(void) {
    quicked_aligner_t aligner;                          // Aligner object
    quicked_status_t status;                            // Return code from QuickEdit functions
    quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters

    params.algo = HIRSCHBERG;                           // Select the algorithm: Hirschberg
    //params.only_score = true;                         // This parameter is ignored when using Hirschberg.

    status = quicked_new(&aligner, &params);            // Initialize the aligner with the given parameters

    const char* pattern = "ACGT";                       // Pattern sequence
    const char* text = "ACTT";                          // Text sequence

    // Align the sequences!
    printf("Aligning '%s' and '%s' using Hirschberg\n", pattern, text);

    status = quicked_align(&aligner, pattern, strlen(pattern), text, strlen(text));
    if (quicked_check_error(status)) {
        fprintf(stderr, "%s", quicked_status_msg(status));
        return 1;
    }

    printf("Score: %d\n", aligner.score);   // Print the score
    printf("CIGAR: %s\n", aligner.cigar);   // Print the CIGAR string

    status = quicked_free(&aligner);        // Free whatever memory the aligner allocated
    if (quicked_check_error(status)) {
        fprintf(stderr, "%s", quicked_status_msg(status));
        return 1;
    }

    return 0;
}