#include "pybind11/pybind11.h"

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
