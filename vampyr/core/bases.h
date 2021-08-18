/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <MRCPP/core/InterpolatingBasis.h>
#include <MRCPP/core/LegendreBasis.h>
#include <MRCPP/core/ScalingBasis.h>

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void bases(py::module &m) {
    py::class_<ScalingBasis>(m, "ScalingBasis")
        .def(py::init<int, int>(), "order"_a, "type"_a)
        .def("getScalingType", &ScalingBasis::getScalingType)
        .def("getScalingOrder", &ScalingBasis::getScalingOrder)
        .def("getQuadratureOrder", &ScalingBasis::getQuadratureOrder)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__repr__", [](const ScalingBasis &basis) {
            std::ostringstream os;
            os << basis;
            return os.str();
        });

    py::class_<InterpolatingBasis, ScalingBasis>(m, "InterpolatingBasis")
        .def(py::init<int>(), "order"_a);

    py::class_<LegendreBasis, ScalingBasis>(m, "LegendreBasis")
        .def(py::init<int>(), "order"_a);
}

} // namespace vampyr
