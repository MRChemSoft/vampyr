/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

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
//#include "functions/AnalyticFunction.h"
#include "functions/BoysFunction.h"
#include "functions/LegendrePoly.h"


using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void bases(py::module &);
void project(py::module &);
void math(py::module &);
void world(py::module &);
void trees(py::module &);
void applys(py::module &);
void derivative(py::module &);
void representable_functions(py::module &);

PYBIND11_MODULE(vampyr3d, m) {

    const auto D = 3;

    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    bases(m);
    project(m);
    math(m);
    world(m);
    trees(m);
    applys(m);
    derivative(m);
    representable_functions(m);





    // m.def("map",
    //       py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, RepresentableFunction<D> &>(&map<D>),
    //       "prec"_a,
    //       "out"_a,
    //       "inp"_a,
    //       "Func"_a);


    m.def(
        "divergence",
        py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTreeVector<D> &>(&mrcpp::divergence<D>));

    m.def("gradient", py::overload_cast<DerivativeOperator<D> &, FunctionTree<D> &>(&gradient<D>));

    py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
    convop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    py::class_<IdentityConvolution<D>> identityconvop(m, "IdentityConvolution");
    identityconvop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    py::class_<PHOperator<D>> phoperator(m, "PHOperator");
    phoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());

    py::class_<BSOperator<D>> bsoperator(m, "BSOperator");
    bsoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());

    py::class_<OperatorState<D>> operatorstate(m, "OperatorState");
    operatorstate.def(py::init<MWNode<D> &, double*>());

    py::class_<GreensKernel> greensKernel(m, "GreensKernel");

    py::class_<MWOperator> mwoperator(m, "MWOperator");
    mwoperator.def(py::init<MultiResolutionAnalysis<2>>());
    
    py::class_<DerivativeKernel>(m, "DerivativeKernel", greensKernel)
        .def(py::init<double>(),
             "eps"_a,
             "GreensKernel: DerivativeKernel");

    py::class_<HelmholtzKernel>(m, "HelmholtzKernel", greensKernel)
        .def(py::init<double, double, double, double>(),
             "m"_a,
             "eps"_a,
             "r_min"_a,
             "r_max"_a,
             "GreensKernel: HelmholtzKernel");    

    py::class_<IdentityKernel>(m, "IdentityKernel", greensKernel)
        .def(py::init<double>(),
             "eps"_a,
             "GreensKernel: IdentityKernel");

    py::class_<PoissonOperator>(m, "PoissonOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, int>(),
             "MRA"_a,
             "precision"_a,
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");

    py::class_<HelmholtzOperator>(m, "HelmholtzOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "MRA"_a,
             "mu"_a,
             "precision"_a,
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");

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
