#pragma once

#include <MRCPP/functions/function_utils.h>
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
        .def("pow", py::overload_cast<int>(&Gaussian<D>::getPower, py::const_), "dim"_a)
        .def("exp", py::overload_cast<int>(&Gaussian<D>::getExp, py::const_), "dim"_a=0)
        .def("pos", &Gaussian<D>::getPos)
        .def("coef", &Gaussian<D>::getCoef)
        .def("periodify", &Gaussian<D>::periodify, "period"_a, "std_dev"_a=4.0)
        .def("overlap", &Gaussian<D>::calcOverlap, "inp"_a)
        .def("__str__", [](const Gaussian<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        });

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
        .def("differentiate", [](const GaussFunc<D> &gauss, int dir) { return gauss.differentiate(dir).asGaussExp(); }, "dir"_a)
        .def("squaredNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);

    // GaussExp class
    py::class_<GaussExp<D>, RepresentableFunction<D>>(m, "GaussExp")
        .def(py::init())
        .def("size", py::overload_cast<>(&GaussExp<D>::size, py::const_))
        .def("func", py::overload_cast<int>(&GaussExp<D>::getFunc), "term"_a, py::return_value_policy::reference_internal)
        .def("append", py::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append))
        .def("periodify", &GaussExp<D>::periodify, "period"_a, "std_dev"_a=4.0)
        .def("differentiate", &GaussExp<D>::differentiate, "dir"_a)
        .def("squaredNorm", &GaussExp<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("__str__", [](const GaussExp<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        });

}
} // namespace vampyr
