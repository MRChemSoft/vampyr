#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/functions/Gaussian.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "AnalyticFunction.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "PyRepresentableFunction.h"

namespace vampyr {
void legendre_poly(pybind11::module &, pybind11::class_<mrcpp::RepresentableFunction<1>, PyRepresentableFunction<1>> &);
void polynomial(pybind11::module &, pybind11::class_<mrcpp::RepresentableFunction<1>, PyRepresentableFunction<1>> &);
void boys_function(pybind11::module &, pybind11::class_<mrcpp::RepresentableFunction<1>, PyRepresentableFunction<1>> &);

template <int D> void functions(pybind11::module &m) {
    using namespace mrcpp;

    pybind11::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, "RepresentableFunction");
    repfunc.def(pybind11::init<>());
    analytic_function(m, repfunc);
    gauss_exp(m, repfunc);

    pybind11::class_<Gaussian<D>> gaussian(m, "Gaussian", repfunc);
    gauss_func(m, gaussian);
    gauss_poly(m, gaussian);
}
} // namespace vampyr
