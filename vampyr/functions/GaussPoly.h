#pragma once

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>

#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/Gaussian.h>

namespace vampyr {
template <int D> void gauss_poly(pybind11::module &m, pybind11::class_<mrcpp::Gaussian<D>> &gaussian) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    pybind11::class_<GaussPoly<D>>(m, "GaussPoly", gaussian)
        .def(pybind11::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("setPoly", &GaussPoly<D>::setPoly)
        .def("differentiate", &GaussPoly<D>::differentiate)
        .def("calcOverlap", pybind11::overload_cast<GaussPoly<D> &>(&GaussPoly<D>::calcOverlap))
        .def("evalf", pybind11::overload_cast<const Coord<D> &>(&GaussPoly<D>::evalf, pybind11::const_))
        .def("evalf", pybind11::overload_cast<const double, int>(&GaussPoly<D>::evalf, pybind11::const_));
}
} // namespace vampyr
