#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>

#include <MRCPP/functions/GaussExp.h>
#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/RepresentableFunction.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D>
void gauss_exp(pybind11::module &m, pybind11::class_<mrcpp::RepresentableFunction<D>, PyRepresentableFunction<D>> &repfunc) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    pybind11::class_<GaussExp<D>>(m, "GaussExp", repfunc)
        .def(pybind11::init<int, double>(), "nTerms"_a = 0, "Gauss_exp_prec"_a = 1.0e-10)
        .def(pybind11::init<const GaussExp<D> &>())
        .def(pybind11::init<const GaussPoly<D> &>())
        .def("evalf", &GaussExp<D>::evalf)
        .def("mult", pybind11::overload_cast<GaussExp<D> &>(&GaussExp<D>::mult))
        .def("mult", pybind11::overload_cast<GaussFunc<D> &>(&GaussExp<D>::mult))
        .def("mult", pybind11::overload_cast<GaussPoly<D> &>(&GaussExp<D>::mult))
        .def("mult", pybind11::overload_cast<double>(&GaussExp<D>::mult))
        .def("multInPlace", pybind11::overload_cast<double>(&GaussExp<D>::multInPlace))
        .def("differentiate", &GaussExp<D>::differentiate)
        .def("getSquareNorm", &GaussExp<D>::getSquareNorm)
        .def("normalize", &GaussExp<D>::normalize)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("add", pybind11::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("add", pybind11::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("append", pybind11::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append))
        .def("append", pybind11::overload_cast<const GaussExp<D> &>(&GaussExp<D>::append));
}
} // namespace vampyr
