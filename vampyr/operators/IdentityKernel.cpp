#include <pybind11/pybind11.h>

#include <MRCPP/operators/GreensKernel.h>
#include <MRCPP/operators/IdentityKernel.h>

namespace vampyr {
using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

void identity_kernel(py::module &m) {
    //IdentityKernel
    py::class_<IdentityKernel, GreensKernel>(m, "IdentityKernel")
        .def(py::init<double>(), "eps"_a, "GreensKernel: IdentityKernel");
}
} // namespace vampyr
