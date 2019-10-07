#include "pybind11/pybind11.h"

#include "operators/IdentityKernel.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void identity_kernel(py::module &m, auto &greenskernel) {

    //IdentityKernel
    py::class_<IdentityKernel>(m, "IdentityKernel", greenskernel)
        .def(py::init<double>(),
             "eps"_a,
             "GreensKernel: IdentityKernel");


}
} // namespace vampyr