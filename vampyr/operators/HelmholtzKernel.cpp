#include <array>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "operators/GreensKernel.h"
#include "operators/HelmholtzKernel.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {


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
