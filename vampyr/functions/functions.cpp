#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "PyRepresentableFunction.h"

#include "functions/Gaussian.h"
#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

    void analytic_function(py::module &, py::class_<RepresentableFunction<3>, PyRepresentableFunction<3>> &);
    void boys_function(py::module &, py::class_<RepresentableFunction<3>, PyRepresentableFunction<3>> &);
    void gauss_exp(py::module &, py::class_<RepresentableFunction<3>, PyRepresentableFunction<3>> &);
    void gauss_func(py::module &, py::class_<Gaussian<3>> &);
    void gauss_poly(py::module &, py::class_<Gaussian<3>> &);
    void legendre_poly(py::module &, py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> &);
    void polynomial(py::module &, py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> &);

void functions(py::module &m) {
    const auto D = 3;


    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, "RepresentableFunction");
    repfunc.def(py::init<>());
    // We need this for the polynomial
    py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> repfunc1d(m, "RepresentableFunction1D");
    repfunc.def(py::init<>());

    py::class_<Gaussian<D>> gaussian(m, "Gaussian", repfunc);

    analytic_function(m, repfunc);
    boys_function(m, repfunc);
    gauss_exp(m, repfunc);
    gauss_func(m, gaussian);
    gauss_poly(m, gaussian);
    legendre_poly(m, repfunc1d);
    polynomial(m, repfunc1d);


}
} // namespace vampyr
