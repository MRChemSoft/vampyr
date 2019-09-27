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
#include "operators/IdentityConvolution.h"
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
void operators(py::module &m) {
    const auto D = 3;


 	  m.def(
        "divergence",
        py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTreeVector<D> &>(&mrcpp::divergence<D>));

    m.def("gradient", py::overload_cast<DerivativeOperator<D> &, FunctionTree<D> &>(&gradient<D>));
 
   //ABGVOperator
    py::class_<ABGVOperator<D>> abgvoperator(m, "ABGVOperator");
    abgvoperator.def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "Derivative Operator: The ABGVOperator for differentiation");

    //BSOperator
 	  py::class_<BSOperator<D>> bsoperator(m, "BSOperator");
    bsoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());

    //ConvolutionOperator
	  py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
    convop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    //DerivativeKernel
    py::class_<DerivativeKernel> derivativekernel(m, "DerivativeKernel");
    derivativekernel.def(py::init<double>(),
             "eps"_a,
             "GreensKernel: DerivativeKernel");

    //DerivativeOperator
    py::class_<DerivativeOperator<D>> derivativeoperator(m, "Derivative Operator");
    derivativeoperator.def(py::init<MultiResolutionAnalysis<D>>());

    //GreensKernel
	  py::class_<GreensKernel> greensKernel(m, "GreensKernel");

    //HelmholtzKernel
    py::class_<HelmholtzKernel>(m, "HelmholtzKernel", greensKernel)
        .def(py::init<double, double, double, double>(),
             "m"_a,
             "eps"_a,
             "r_min"_a,
             "r_max"_a,
             "GreensKernel: HelmholtzKernel");    

    //HelmholtzOperator
    py::class_<HelmholtzOperator, ConvolutionOperator<D>, MWOperator>(m, "HelmholtzOperator", py::multiple_inheritance())
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "MRA"_a,
             "mu"_a,
             "precision"_a,
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");

    //IdentityConvolution
	  py::class_<IdentityConvolution<D>> identityconvop(m, "IdentityConvolution");
    identityconvop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    //IdentityKernel
    py::class_<IdentityKernel>(m, "IdentityKernel", greensKernel)
        .def(py::init<double>(),
             "eps"_a,
             "GreensKernel: IdentityKernel");

    //MWOperator
    py::class_<MWOperator> mwoperator(m, "MWOperator");
    mwoperator.def(py::init<MultiResolutionAnalysis<2>>())
        .def("size", &MWOperator::size, "Size")
        .def("push_back", &MWOperator::push_back, "Push")
        .def("clear", &MWOperator::clear, "Clear")
        .def("getMaxBandWidth", &MWOperator::getMaxBandWidth, "Get")
        .def("calcBandWidths", &MWOperator::calcBandWidths, "Calc")
        .def("clearBandWidths", &MWOperator::clearBandWidths, "Clear");

    //OperatorState
    py::class_<OperatorState<D>> operatorstate(m, "OperatorState");
    operatorstate.def(py::init<MWNode<D> &, double*>());

    //OperatorStatistics
    py::class_<OperatorStatistics<D>> operatorstatistics(m, "OperatorStatistics");
    operatorstatistics.def(py::init<>())
        .def("flushNodeCounters", &OperatorStatistics<D>::flushNodeCounters, "flushNodeCounters")
        .def("incrementFNodeCounters", &OperatorStatistics<D>::incrementFNodeCounters, "incrementFNodeCounters")
        .def("incrementGNodeCounters", &OperatorStatistics<D>::incrementGNodeCounters, "incrementGNodeCounters");

    //PHOperator
    py::class_<PHOperator<D>> phoperator(m, "PHOperator");
    phoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());

    //PoissonOperator
    py::class_<PoissonOperator, ConvolutionOperator<D>, MWOperator>(m, "PoissonOperator", py::multiple_inheritance())
        .def(py::init<MultiResolutionAnalysis<D> &, double>(),
             "MRA"_a,
             "precision"_a,
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");


    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, const RepresentableFunction<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTreeVector<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_grid<D>));
    m.def("copy_func", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_func<D>));

    m.def("clear_grid", py::overload_cast<FunctionTree<D> &>(&clear_grid<D>));
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, int>(&refine_grid<D>), "out"_a, "scales"_a);

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, double>(&refine_grid<D>), "out"_a, "prec"_a);
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&refine_grid<D>), "out"_a, "inp"_a);

}
} // namespace vampyr
