/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <MRCPP/core/InterpolatingBasis.h>
#include <MRCPP/core/LegendreBasis.h>
#include <MRCPP/core/ScalingBasis.h>

#include "functions/PyRepresentableFunction.h"
#include "functions/PyScalingFunction.h"
#include "functions/PyWaveletFunction.h"

namespace vampyr {

void bases(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<Polynomial, PyRepresentableFunction<1, Polynomial>, RepresentableFunction<1, double>>(m, "Polynomial");

    py::class_<PyScalingFunction, Polynomial>(m, "ScalingFunction");

    py::class_<PyWaveletFunction, RepresentableFunction<1, double>>(m, "WaveletFunction");

    py::class_<ScalingBasis>(m,
                             "ScalingBasis",
                             R"mydelimiter(
        Filler text
    )mydelimiter")
        .def(py::init<int, int>(), "order"_a, "type"_a)
        .def("scalingType", &ScalingBasis::getScalingType)
        .def("scalingOrder", &ScalingBasis::getScalingOrder)
        .def("quadratureOrder", &ScalingBasis::getQuadratureOrder)
        .def(
            "scaling",
            [](const ScalingBasis &basis, int i, int l, int n) { return PyScalingFunction(basis, i, l, n); },
            "i"_a,
            "l"_a = 0,
            "n"_a = 0)
        .def(
            "wavelet",
            [](const ScalingBasis &basis, int i, int l, int n) { return PyWaveletFunction(basis, i, l, n); },
            "i"_a,
            "l"_a = 0,
            "n"_a = 0)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__str__", [](const ScalingBasis &basis) {
            std::ostringstream os;
            os << basis;
            return os.str();
        });

    py::class_<InterpolatingBasis, ScalingBasis>(m, "InterpolatingBasis").def(py::init<int>(), "order"_a);

    py::class_<LegendreBasis, ScalingBasis>(m, "LegendreBasis").def(py::init<int>(), "order"_a);
}

} // namespace vampyr
