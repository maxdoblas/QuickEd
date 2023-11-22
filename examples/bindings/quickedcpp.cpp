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

#include <quicked.hpp>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(void) {
    quicked::QuickedAligner aligner;                    // Aligner object, with sensible default parameters
    // Without any extra configuration, the aligner will use the Quicked algorithm

    string pattern = "ACGT";                            // Pattern sequence
    string text = "ACTT";                               // Text sequence

    // Align the sequences!
    cout << "Aligning " << pattern << " and " << text << " using Quicked" << endl;
    aligner.align(&pattern, pattern.length(), &text, text.length());

    cout << "Score: " << aligner.get_score() << endl;   // Print the score
    cout << "Cigar: " << aligner.get_cigar() << endl;   // Print the CIGAR string

    return 0;
}