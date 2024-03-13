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

#include <pybind11/pybind11.h>
#include "quicked.hpp"

namespace py = pybind11;
namespace quicked {
    PYBIND11_MODULE(pyquicked, m) {
        m.doc() = "QuickEd Library";

    py::class_<QuickedAligner>(m, "QuickedAligner")
            .def(py::init<>())
            .def("align", &QuickedAligner::align)
            .def("setAlgorithm", &QuickedAligner::setAlgorithm)
            .def("setOnlyScore", &QuickedAligner::setOnlyScore)
            .def("setBandwidth", &QuickedAligner::setBandwidth)
            .def("setWindowSize", &QuickedAligner::setWindowSize)
            .def("setOverlapSize", &QuickedAligner::setOverlapSize)
            .def("setForceScalar", &QuickedAligner::setForceScalar)
            .def("setHEWThreshold", &QuickedAligner::setHEWThreshold)
            .def("setHEWPercentage", &QuickedAligner::setHEWPercentage)
            .def("getScore", &QuickedAligner::getScore)
            .def("getCigar", &QuickedAligner::getCigar);

        py::enum_<quicked_algo_t>(m, "QuickedAlgo")
            .value("QUICKED", QUICKED)
            .value("WINDOWED", WINDOWED)
            .value("BANDED", BANDED)
            .value("HIRSCHBERG", HIRSCHBERG)
            .export_values();

        py::enum_<quicked_status_t>(m, "QuickedStatus")
            .value("QUICKED_OK", QUICKED_OK)
            .value("QUICKED_ERROR", QUICKED_ERROR)
            .value("QUICKED_FAIL_NON_CONVERGENCE", QUICKED_FAIL_NON_CONVERGENCE)
            .value("QUICKED_UNKNOWN_ALGO", QUICKED_UNKNOWN_ALGO)
            .value("QUICKED_EMPTY_SEQUENCE", QUICKED_EMPTY_SEQUENCE)
            .value("QUICKED_UNIMPLEMENTED", QUICKED_UNIMPLEMENTED)
            .value("QUICKED_WIP", QUICKED_WIP)
            .export_values();

        py::register_exception<QuickedException>(m, "QuickedException");

    }
}