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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char *argv[]) { // Usage: ./quicked_harness <text> <pattern> <expected score>

    quicked_aligner_t aligner;
    quicked_status_t status;
    quicked_params_t params = quicked_default_params();

    status = quicked_new(&aligner, &params);
    if (quicked_check_error(status)) {
        fprintf(stderr, "%s", quicked_status_msg(status));
        exit(EXIT_FAILURE);
    }

    status = quicked_align(&aligner, argv[1], strlen(argv[1]), argv[2], strlen(argv[2]));
    if (quicked_check_error(status)) {
        fprintf(stderr, "%s", quicked_status_msg(status));
        exit(EXIT_FAILURE);
    }

    int score = aligner.score;
    printf("Got score: %d\n", score);

    status = quicked_free(&aligner);
    if (quicked_check_error(status)) {
        fprintf(stderr, "%s", quicked_status_msg(status));
        exit(EXIT_FAILURE);
    }

    if (argc == 4) { // If expected score is provided, check if it matches
        printf("Expected score: %d", atoi(argv[3]));
        if (score != atoi(argv[3])) {
            printf("<FAIL>\n");
            exit(EXIT_FAILURE);
        }
    }

    return 0;
}