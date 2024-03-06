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

namespace quicked {

    QuickedAligner::QuickedAligner()
    {
        this->params = quicked_default_params();
        quicked_new(&this->aligner, &this->params);
    }

    QuickedAligner::~QuickedAligner()
    {
        quicked_free(&this->aligner);
    }

    quicked_status_t QuickedAligner::align(std::string *pattern, std::string *text)
    {
        return quicked_align(&this->aligner, pattern->c_str(), pattern->length(), text->c_str(), text->length());
    }

    void QuickedAligner::setHEWThreshold(unsigned int hew_threshold) {
        for (int i = 0; i < QUICKED_WINDOW_STAGES; i++) this->aligner.params->hew_threshold[i] = hew_threshold;
    };
    void QuickedAligner::setHEWPercentage(unsigned int hew_percentage) {
        for (int i = 0; i < QUICKED_WINDOW_STAGES; i++) this->aligner.params->hew_percentage[i] = hew_percentage;
    };
}