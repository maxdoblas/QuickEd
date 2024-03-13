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

#include <quicked.hpp>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(void) {
    string pattern = "ACGT"; // Pattern sequence
    string text = "ACTT";    // Text sequence
    int score = -1;          // Alignment score
    string cigar;            // CIGAR string

    cout << "Aligning " << pattern << " and " << text << " using Quicked" << endl;

    try {
        quicked::QuickedAligner aligner; // Aligner object, with sensible default parameters
        // Without any extra configuration, the aligner will use the Quicked algorithm

        aligner.align(&pattern, &text);  // Align the sequences!

        aligner.getScore();              // Get the score
        aligner.getCigar();              // Get the CIGAR string
    } catch (quicked::QuickedException &e) {
        cerr << e.what() << endl;
        return 1;
    }

    cout << "Score: " << score << endl;  // Print the score
    cout << "Cigar: " << cigar << endl;  // Print the CIGAR string

    return 0;
}