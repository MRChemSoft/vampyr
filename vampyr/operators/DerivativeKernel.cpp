#include <array>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "operators/DerivativeKernel.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void derivative_kernel(py::module &m) {
    const auto D = 3;

    //DerivativeKernel
    py::class_<DerivativeKernel> derivativekernel(m, "DerivativeKernel");
    derivativekernel.def(py::init<double>(),
             "eps"_a,
             "GreensKernel: DerivativeKernel");


}
} // namespace vampyr
