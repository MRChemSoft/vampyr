#pragma once

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/functions/Gaussian.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "AnalyticFunction.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "PyRepresentableFunction.h"

namespace vampyr {
void polynomial(pybind11::module &);
void legendre_poly(pybind11::module &);
void boys_function(pybind11::module &);

template <int D> void functions(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>>(m, "RepresentableFunction")
        .def(py::init<const std::vector<double> &, const std::vector<double> &>(), "a"_a, "b"_a)
        .def("evalf", &RepresentableFunction<D>::evalf);

    analytic_function<D>(m);
    gauss_exp<D>(m);

    py::class_<Gaussian<D>, RepresentableFunction<D>>(m, "Gaussian");
    // FIXME These cannot currently be bound due to confusion with inheritance hierarchies.
    //    .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>(), "a"_a, "c"_a, "r"_a, "p"_a)
    //.def(py::init<const std ::array<double, D> &, double, const Coord<D> &, const std::array<int, D> &>(),
    //     "a"_a,
    //     "c"_a,
    //     "r"_a,
    //     "p"_a)
    //.def("evalf", py::overload_cast<const Coord<D>&>(&Gaussian<D>::evalf, py::const_))
    // Bind the other crap
    //.def("__str__", [](const Gaussian<D> &obj) {
    //    std::ostringstream os;
    //    os << obj;
    //    return os.str();
    //});
    gauss_func<D>(m);
    gauss_poly<D>(m);
}
} // namespace vampyr
