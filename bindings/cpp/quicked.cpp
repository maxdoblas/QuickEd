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

#include "quicked.hpp"

quicked::QuickedAligner::QuickedAligner()
{
    quicked_new(&(this->aligner), quicked_default_params());
}

quicked::QuickedAligner::QuickedAligner(quicked_params_t params)
{
    quicked_new(&(this->aligner), params);
}

quicked::QuickedAligner::~QuickedAligner()
{
    quicked_free(&(this->aligner));
}

quicked_status_t quicked::QuickedAligner::align(
    std::string *pattern, const int pattern_len,
    std::string *text, const int text_len)
{
    return quicked_align(&(this->aligner), pattern->c_str(), pattern_len, text->c_str(), text_len);
}
