#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/GaussExp.h>
#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D> void gauss_exp(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<GaussExp<D>, RepresentableFunction<D>>(m, "GaussExp")
        .def(py::init<int, double>(), "nTerms"_a = 0, "Gauss_exp_prec"_a = 1.0e-10)
        .def(py::init<const GaussExp<D> &>())
        .def(py::init<const GaussPoly<D> &>())
        .def("mult", py::overload_cast<GaussExp<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<GaussFunc<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<GaussPoly<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<double>(&GaussExp<D>::mult))
        .def("multInPlace", py::overload_cast<double>(&GaussExp<D>::multInPlace))
        .def("differentiate", &GaussExp<D>::differentiate)
        .def("getSquareNorm", &GaussExp<D>::getSquareNorm)
        .def("normalize", &GaussExp<D>::normalize)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("add", py::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("add", py::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("append", py::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append))
        .def("append", py::overload_cast<const GaussExp<D> &>(&GaussExp<D>::append));
}
} // namespace vampyr
