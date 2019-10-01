#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "trees/MWTree.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/OperatorTree.h"
#include "trees/NodeIndex.h"
#include "trees/MWNode.h"
#include "trees/BandWidth.h"
#include "trees/FunctionNode.h"
#include "trees/OperatorNode.h"
#include "trees/SerialTree.h"
#include "trees/SerialFunctionTree.h"
#include "trees/SerialOperatorTree.h"
#include "trees/TreeIterator.h"

#include "utils/mpi_utils.h"

#include "PyRepresentableFunction.h"

#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "operators/HelmholtzOperator.h"
#include "operators/PoissonOperator.h"
#include "operators/PHOperator.h"
#include "operators/BSOperator.h"
#include "operators/GreensKernel.h"
#include "operators/DerivativeKernel.h"
#include "operators/IdentityKernel.h"
#include "operators/MWOperator.h"
#include "operators/HelmholtzKernel.h"
#include "operators/OperatorState.h"
#include "operators/OperatorStatistics.h"

#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

#include "functions/GaussFunc.h"
#include "functions/GaussExp.h"
#include "functions/GaussPoly.h"
#include "functions/Gaussian.h"
#include "functions/Polynomial.h"
#include "functions/RepresentableFunction.h"
#include "functions/AnalyticFunction.h"
#include "functions/BoysFunction.h"
#include "functions/LegendrePoly.h"


using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void functions(py::module &m) {
    const auto D = 3;


    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, "RepresentableFunction");
    repfunc.def(py::init<>());
    // We need this for the polynomial
    py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> repfunc1d(m, "RepresentableFunction1D");
    repfunc.def(py::init<>());

    py::class_<AnalyticFunction<D>>(m, "AnalyticFunction", repfunc)
        .def(py::init<std::function<double(const Coord<D>&r)>, double *, double *>(), "r"_a = Coord<D>{}, "a"_a = nullptr, "b"_a = nullptr)
        .def("evalf", py::overload_cast<const Coord<D> &>(&AnalyticFunction<D>::evalf, py::const_), 
         "This function does not work properly since the input are raw pointers, for this function to work properly add input arrays in MRCPP");

    py::class_<BoysFunction>(m, "BoysFunction", repfunc)
        .def(py::init<int, double>(), "nTerms"_a = 0, "Boys_func_prec"_a = 1.0e-10)
        .def("evalf", py::overload_cast<const Coord<1> &>(&BoysFunction::evalf, py::const_));

    py::class_<Gaussian<D>> gaussian(m, "Gaussian", repfunc);

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

    py::class_<GaussFunc<D>>(m, "GaussFunc", gaussian)
        .def(py::init<double, double, Coord<D> &, std::array<int, D> &>(), "alpha"_a, "beta"_a, "pos"_a = Coord<D>{}, "power"_a = std::array<int, D>{})
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussFunc<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<double, int>(&GaussFunc<D>::evalf, py::const_))
        .def("calcOverlap", py::overload_cast<GaussFunc<D> &>(&GaussFunc<D>::calcOverlap))
        .def("differentiate", &GaussFunc<D>::differentiate)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy)
        .def("getPos", &GaussFunc<D>::getPos)
        .def("getExp", &GaussFunc<D>::getExp)
        .def("getCoef", &GaussFunc<D>::getCoef)
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);

    py::class_<GaussPoly<D>>(m, "GaussPoly", gaussian)
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("setPoly", &GaussPoly<D>::setPoly)
        .def("differentiate", &GaussPoly<D>::differentiate)
        .def("calcOverlap", py::overload_cast<GaussPoly<D> &>(&GaussPoly<D>::calcOverlap))
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussPoly<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<const double, int>(&GaussPoly<D>::evalf, py::const_));

    py::class_<LegendrePoly>(m, "LegendrePoly", repfunc1d)
        .def(py::init<int, double, double>(), "k"_a, "n"_a = 1.0, "l"_a = 0.0)
        .def("firstDerivative", &LegendrePoly::firstDerivative)
        .def("secondDerivative", &LegendrePoly::secondDerivative);

    py::class_<Polynomial>(m, "Polynomial", repfunc1d)
        .def(py::init<int, const double *, const double *>(), "order"_a, "a"_a = nullptr, "b"_a = nullptr)
        .def("evalf", py::overload_cast<double>(&Polynomial::evalf, py::const_))
        .def("setCoefs", &Polynomial::setCoefs)
        .def("normalize", &Polynomial::normalize)
        .def("getCoefs", py::overload_cast<>(&Polynomial::getCoefs))
        .def("size", &Polynomial::size);

}
} // namespace vampyr
