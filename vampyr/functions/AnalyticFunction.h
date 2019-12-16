#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/AnalyticFunction.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D>
void analytic_function(pybind11::module &m,
                       pybind11::class_<mrcpp::RepresentableFunction<D>, PyRepresentableFunction<D>> &repfunc) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // The double * wont work properly in vampyr since raw pointers does not exist in python.
    // This should probably be changed to an std::vector in mrcpp.
    py::class_<AnalyticFunction<D>>(m, "AnalyticFunction", repfunc)
        .def(py::init<std::function<double(const Coord<D> &r)>, double *, double *>(),
             "f"_a,
             "a"_a = nullptr,
             "b"_a = nullptr)
        .def("evalf",
             py::overload_cast<const Coord<D> &>(&AnalyticFunction<D>::evalf, py::const_),
             "This function does not work properly since the input are raw pointers, for this function to work "
             "properly add input arrays in MRCPP");
}
} // namespace vampyr
