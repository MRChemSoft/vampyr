#include <pybind11/pybind11.h>

#include <MRCPP/operators/GreensKernel.h>
#include <MRCPP/operators/HelmholtzKernel.h>

namespace vampyr {
using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

void helmholtz_kernel(py::module &m, py::class_<GreensKernel> &greenskernel) {
    //HelmholtzKernel
    py::class_<HelmholtzKernel>(m, "HelmholtzKernel", greenskernel)
        .def(py::init<double, double, double, double>(),
             "m"_a,
             "eps"_a,
             "r_min"_a,
             "r_max"_a,
             "GreensKernel: HelmholtzKernel");
}
} // namespace vampyr
