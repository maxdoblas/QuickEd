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

#include <pybind11/pybind11.h>
#include "quicked.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyquicked, m) {
    m.doc() = "QuickEd Library";

py::class_<quicked::QuickedAligner>(m, "QuickedAligner")
        .def(py::init<>())
        .def("align", &quicked::QuickedAligner::align)
        .def("setAlgorithm", &quicked::QuickedAligner::setAlgorithm)
        .def("setOnlyScore", &quicked::QuickedAligner::setOnlyScore)
        .def("setBandwidth", &quicked::QuickedAligner::setBandwidth)
        .def("setWindowSize", &quicked::QuickedAligner::setWindowSize)
        .def("setOverlapSize", &quicked::QuickedAligner::setOverlapSize)
        .def("setForceScalar", &quicked::QuickedAligner::setForceScalar)
        .def("getScore", &quicked::QuickedAligner::getScore)
        .def("getCigar", &quicked::QuickedAligner::getCigar);

    py::enum_<quicked::quicked_algo_t>(m, "QuickedAlgo")
        .value("QUICKED", quicked::QUICKED)
        .value("WINDOWED", quicked::WINDOWED)
        .value("BANDED", quicked::BANDED)
        .value("HIRSCHBERG", quicked::HIRSCHBERG)
        .export_values();

    py::enum_<quicked::quicked_status_t>(m, "QuickedStatus")
        .value("QUICKED_OK", quicked::QUICKED_OK)
        .value("QUICKED_ERROR", quicked::QUICKED_ERROR)
        .value("QUICKED_UNKNOWN_ALGO", quicked::QUICKED_UNKNOWN_ALGO)
        .value("QUICKED_UNIMPLEMENTED", quicked::QUICKED_UNIMPLEMENTED)
        .value("QUICKED_WIP", quicked::QUICKED_WIP)
        .export_values();

}