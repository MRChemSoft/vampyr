#include "pybind11/pybind11.h"

#include "PyRepresentableFunction.h"

#include "functions/GaussPoly.h"
#include "functions/GaussExp.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void gauss_exp(py::module &m, py::class_<RepresentableFunction<3>, PyRepresentableFunction<3>> &repfunc) {
    
    const auto D = 3;


py::class_<GaussExp<D>>(m, "GaussExp", repfunc)
        .def(py::init<int, double>(), "nTerms"_a = 0, "Gauss_exp_prec"_a = 1.0e-10)
        .def(py::init<const GaussExp<D> &>())
        .def(py::init<const GaussPoly<D> &>())
        .def("evalf", &GaussExp<D>::evalf)
        .def("mult", py::overload_cast<GaussExp<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<GaussFunc<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<GaussPoly<D> &>(&GaussExp<D>::mult))
        .def("mult", py::overload_cast<double>(&GaussExp<D>::mult))
        .def("multInPlace", py::overload_cast<double>(&GaussExp<D>::multInPlace))
        .def("differentiate", &GaussExp<D>::differentiate)
        .def("getSquareNorm", &GaussExp<D>::getSquareNorm)
        .def("normalize", &GaussExp<D>::normalize)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("add", py::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("add", py::overload_cast<Gaussian<D> &>(&GaussExp<D>::add))
        .def("append", py::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append))
        .def("append", py::overload_cast<const GaussExp<D> &>(&GaussExp<D>::append));


}
} // namespace vampyr
