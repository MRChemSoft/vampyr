#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "functions/Gaussian.h"
#include "functions/GaussPoly.h"
#include "functions/GaussFunc.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void gauss_func(py::module &m, py::class_<Gaussian<3>> &gaussian) {
    const auto D = 3;


py::class_<GaussFunc<D>>(m, "GaussFunc", gaussian)
        .def(py::init<double, double, Coord<D> &, std::array<int, D> &>(), "alpha"_a, "beta"_a, "pos"_a = Coord<D>{}, "power"_a = std::array<int, D>{})
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussFunc<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<double, int>(&GaussFunc<D>::evalf, py::const_))
        .def("calcOverlap", py::overload_cast<GaussFunc<D> &>(&GaussFunc<D>::calcOverlap))
        .def("differentiate", &GaussFunc<D>::differentiate)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy)
        .def("getPos", &GaussFunc<D>::getPos)
        .def("getExp", &GaussFunc<D>::getExp)
        .def("getCoef", &GaussFunc<D>::getCoef)
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);


}
} // namespace vampyr
