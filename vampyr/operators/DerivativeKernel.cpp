#include <pybind11/pybind11.h>

#include <MRCPP/operators/DerivativeKernel.h>
#include <MRCPP/operators/GreensKernel.h>

namespace vampyr {
using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

void derivative_kernel(py::module &m) {
    //DerivativeKernel
    py::class_<DerivativeKernel, GreensKernel>(m, "DerivativeKernel")
        .def(py::init<double>(), "eps"_a, "GreensKernel: DerivativeKernel");
}
} // namespace vampyr
