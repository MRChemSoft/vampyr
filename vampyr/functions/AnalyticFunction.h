#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>

#include <MRCPP/functions/AnalyticFunction.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D>
void analytic_function(pybind11::module &m,
                       pybind11::class_<mrcpp::RepresentableFunction<D>, PyRepresentableFunction<D>> &repfunc) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    // The double * wont work properly in vampyr since raw pointers does not exist in python.
    // This should probably be changed to an std::vector in mrcpp.
    pybind11::class_<AnalyticFunction<D>>(m, "AnalyticFunction", repfunc)
        .def(pybind11::init<std::function<double(const Coord<D> &r)>, double *, double *>(),
             "f"_a,
             "a"_a = nullptr,
             "b"_a = nullptr)
        .def("evalf",
             pybind11::overload_cast<const Coord<D> &>(&AnalyticFunction<D>::evalf, pybind11::const_),
             "This function does not work properly since the input are raw pointers, for this function to work "
             "properly add input arrays in MRCPP");
}
} // namespace vampyr
