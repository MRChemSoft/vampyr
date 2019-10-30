#include <array>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "operators/GreensKernel.h"
#include "operators/IdentityKernel.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void identity_kernel(py::module &m, py::class_<GreensKernel> &greenskernel) {

    //IdentityKernel
    py::class_<IdentityKernel>(m, "IdentityKernel", greenskernel)
        .def(py::init<double>(),
             "eps"_a,
             "GreensKernel: IdentityKernel");


}
} // namespace vampyr
