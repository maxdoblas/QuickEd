#!/usr/bin/env python3

#                             The MIT License
#
#  This file is part of QuickEdit library.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

from pyquicked import QuickedAligner

aligner = QuickedAligner()  # Aligner object, with sensible default parameters
# Without any extra configuration, the aligner will use the Quicked algorithm

pattern = "ACGT";           # Pattern sequence
text = "ACTT";              # Text sequence

#Align the sequences!
print(f"Aligning {pattern} and {text} using Quicked")
aligner.align(pattern, text)

print(f"Score: {aligner.getScore()}") # Print the score
print(f"Cigar: {aligner.getCigar()}")  # Print the CIGAR string
