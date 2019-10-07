#include "pybind11/pybind11.h"

#include "operators/HelmholtzKernel.h"


using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {


void helmholtz_kernel(py::module &m, auto &greenskernel) {

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