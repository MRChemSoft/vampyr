#pragma once

#include <MRCPP/functions/GaussExp.h>
#include <MRCPP/functions/GaussFunc.h>
#include <MRCPP/functions/GaussPoly.h>

#include "PyGaussian.h"

namespace vampyr {

template <int D> void gaussians(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // Gaussian class
    py::class_<Gaussian<D>, PyGaussian<D>, RepresentableFunction<D>>(m, "Gaussian")
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("__str__", [](const Gaussian<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        })
        .def("getPow", py::overload_cast<int>(&Gaussian<D>::getPow, py::const_), "dim"_a)
        .def("getExp", py::overload_cast<int>(&Gaussian<D>::getExp, py::const_), "dim"_a=0)
        .def("getPos", &Gaussian<D>::getPos)
        .def("getCoef", &Gaussian<D>::getCoef);

    // GaussFunc class
    py::class_<GaussFunc<D>, PyGaussian<D, GaussFunc<D>>, Gaussian<D>>(m, "GaussFunc")
        .def(py::init<double,
                      double,
                      const Coord<D> &,
                      const std::array<int, D> &>(),
                      "exp"_a,
                      "coef"_a = 1.0,
                      "pos"_a = Coord<D>{},
                      "pow"_a = std::array<int, D>{})
        .def("calcOverlap", py::overload_cast<GaussFunc<D> &>(&GaussFunc<D>::calcOverlap))
        .def("differentiate", &GaussFunc<D>::differentiate, "dim"_a)
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);

    // GaussPoly class
    py::class_<GaussPoly<D>, PyGaussian<D, GaussPoly<D>>, Gaussian<D>>(m, "GaussPoly");

    // GaussExp class
    py::class_<GaussExp<D>, RepresentableFunction<D>>(m, "GaussExp")
        .def(py::init())
        .def("size", py::overload_cast<>(&GaussExp<D>::size, py::const_))
        .def("append", py::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append))
        .def("getFunc", py::overload_cast<int>(&GaussExp<D>::getFunc), "term"_a, py::return_value_policy::reference_internal)
        .def("differentiate", &GaussExp<D>::differentiate, "dim"_a)
        .def("calcSquareNorm", &GaussExp<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("__str__", [](const GaussExp<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        });
}
} // namespace vampyr
