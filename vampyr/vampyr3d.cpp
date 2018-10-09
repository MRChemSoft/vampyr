/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/MWTree.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "PyAnalyticFunction.h"
#include "PyRepresentableFunction.h"
#include "project.h"

#include "operators/DerivativeOperator.h"
#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/PoissonOperator.h"
#include "operators/HelmholtzOperator.h"

#include "treebuilders/add.h"
#include "treebuilders/multiply.h"
#include "treebuilders/apply.h"
#include "treebuilders/project.h"
#include "treebuilders/grid.h"

#include "functions/GaussFunc.h"
#include "functions/RepresentableFunction.h"
#include "functions/Gaussian.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void bases(py::module &);

PYBIND11_MODULE(vampyr3d, m) {

    const auto D = 3;

    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    bases(m);

    m.def("project", py::overload_cast<double, FunctionTree<D> &,
          std::function<double (std::array<double, D>)>, int>
          (&vampyr::project<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("function"),
          py::arg("maxIter")= -1,
          "Projects an analytic function onto a FunctionTree");

    m.def("project", py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int>(&mrcpp::project<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("GaussFunc"), py::arg("maxIter") = -1);

    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, "RepresentableFunction");
    repfunc
        .def(py::init<>());

    py::class_<Gaussian<D>> gaussian(m, "Gaussian", repfunc);

    py::class_<GaussFunc<D>> (m, "GaussFunc", gaussian)
        .def(py::init<double, double, std::array<double, D> &,
             std::array<int, D> &>())
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);


    py::class_<BoundingBox<D>> (m, "BoundingBox")
        .def(py::init<int, std::array<int, D>, std::array <int, D>>(),
             "Initiates the BoundingBox of the area functions can exist")
        .def("getScale", &BoundingBox<D>::getScale);

    py::class_<MultiResolutionAnalysis<D>> (m, "MultiResolutionAnalysis")
        .def(py::init<BoundingBox<D>, ScalingBasis, int>(),
             py::arg("BoundingBox"), py::arg("Basis Function"),
             py::arg("Scale"))
        .def("getOrder", &MultiResolutionAnalysis<D>::getOrder,
             "Returns the order of the scaling basis")
        .def("getMaxDepth", &MultiResolutionAnalysis<D>::getMaxDepth)
        .def("getMaxScale", &MultiResolutionAnalysis<D>::getMaxScale);

    py::class_<MWTree<D>> mwtree(m, "MWTree");
    mwtree
        .def(py::init<MultiResolutionAnalysis<D>>())
        .def("getSquareNorm", &MWTree<D>::getSquareNorm);

    py::class_<FunctionTree<D>> functree(m, "FunctionTree", mwtree);
    functree
        .def(py::init<MultiResolutionAnalysis<D>>())
        .def("integrate", &FunctionTree<D>::integrate,
             "Integral of the FunctionTree over all space")
        .def("clear", &FunctionTree<D>::clear,
             "Clears the FunctionTree")
        .def("normalize", &FunctionTree<D>::normalize,
             "Rescale the function by its norm, fixed grid")
        .def("rescale", &FunctionTree<D>::rescale,
             "Rescales the function")
        .def("evalf", py::overload_cast<const std::array<double, D> &>
            (&FunctionTree<D>::evalf),
             "Returns the function value at a given point");


    m.def("apply", py::overload_cast<double, FunctionTree<D> &,
          ConvolutionOperator<D> &, FunctionTree<D> &, int>(&apply<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("oper"),
          py::arg("input_tree"), py::arg("maxIter") = -1,
          "Applies a given convolution operator onto a FunctionTree");

    m.def("apply", py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &,
          FunctionTree<D> &, int>(&apply<D>),
          py::arg("out_tree"), py::arg("DerivativeOperator"),
          py::arg("input_tree"), py::arg("dir") = -1,
          "Applies a given derivative operator onto a FunctionTree");

    m.def("add", py::overload_cast<double, FunctionTree<D> &, double,
          FunctionTree<D> &, double, FunctionTree<D> &, int>(&add<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("a"),
          py::arg("tree_a"), py::arg("b"), py::arg("tree_b"),
          py::arg("maxIter") = -1,
          "Addition: output_tree = a*tree_a + b*tree_b");

    m.def("add", py::overload_cast<double, FunctionTree<D> &,
          FunctionTreeVector<D> &, int>(&add<D>),
          py::arg("precision"), py::arg("output_tree"),
          py::arg("FunctionTreeVector"), py::arg("maxIter") = -1,
          "Addition: output_tree = sum_i c_i f_i(x)");

    m.def("multiply", py::overload_cast<double, FunctionTree<D> &, double,
          FunctionTree<D> &, FunctionTree<D> &, int >(&multiply<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("c"),
          py::arg("tree_a"), py::arg("tree_b"), py::arg("maxIter") = -1,
          "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("multiply", py::overload_cast<double, FunctionTree<D> &,
          FunctionTreeVector<D> &, int>(&multiply<D>),
          py::arg("precision"), py::arg("output_tree"), py::arg("inpit_vector"),
          py::arg("maxIter") = -1,
          "Multiplication: output_tree = prod_i c_i f_i(x)");

    m.def("dot", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&dot<D>),
          py::arg("bra"), py::arg("ket"));

    m.def("dot", py::overload_cast<double, FunctionTree<D> &,
          FunctionTreeVector<D> &, FunctionTreeVector<D> &, int>(&dot<D>),
          py::arg("prec"), py::arg("out"), py::arg("inp_a"), py::arg("inp_b"),
          py::arg("maxIter") = -1);

    m.def("power", py::overload_cast<double, FunctionTree<D> &,
          FunctionTree<D> &, double, int>(&power<D>),
          py::arg("prec"), py::arg("out"), py::arg("inp"), py::arg("power"),
          py::arg("maxIter") = -1);

    m.def("square", py::overload_cast<double, FunctionTree<D> &,
          FunctionTree<D> &, int>(&square<D>),
          py::arg("prec"), py::arg("out"), py::arg("inp"),
          py::arg("maxIter") = -1);

    m.def("map", py::overload_cast<double, FunctionTree<D> &,
          FunctionTree<D> &, RepresentableFunction<D> &>(&map<D>),
          py::arg("prec"), py::arg("out"), py::arg("inp"), py::arg("Func"));

    py::class_<DerivativeOperator<D>> deriv(m, "Derivative Operator");
    deriv
        .def(py::init<MultiResolutionAnalysis<D>>());

    py::class_<ABGVOperator<D>> (m, "ABGVOperator", deriv)
        .def(py::init< MultiResolutionAnalysis<D> &, double, double>(),
            "Derivative Operator: The ABGVOperator for differentiation");

    // This does not work outside of jupyter notebooks
    m.def("divergence", py::overload_cast<FunctionTree<D> &,
          DerivativeOperator<D> &, FunctionTreeVector<D> &>
          (&mrcpp::divergence<D>));

    // This does not work outside of jupyter notebooks
    m.def("gradient", py::overload_cast<DerivativeOperator<D> &,
          FunctionTree<D> &>(&gradient<D>));

    py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
    convop
        .def(py::init<MultiResolutionAnalysis<D> &, double>());

    py::class_<PoissonOperator> (m, "PoissonOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, double >(),
             py::arg("MRA"), py::arg("precision"),
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");

    py::class_<HelmholtzOperator> (m, "HelmholtzOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             py::arg("MRA"), py::arg("mu"), py::arg("precision"),
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");

    m.def("build_grid", py::overload_cast<FunctionTree<D> &,
          const RepresentableFunction<D> &, int>(&build_grid<D>),
          py::arg("out"), py::arg("inp"), py::arg("maxIter") = -1);

     m.def("build_grid", py::overload_cast<FunctionTree<D> &,
           FunctionTree<D> &, int>(&build_grid<D>),
           py::arg("out"), py::arg("inp"), py::arg("maxIter") = -1);

    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>
          (&copy_grid<D>));
    m.def("copy_func", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>
          (&copy_func<D>));

    m.def("clear_grid", py::overload_cast<FunctionTree<D> &>(&clear_grid<D>));
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, int>
          (&refine_grid<D>), py::arg("out"), py::arg("scales"));

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, double>
          (&refine_grid<D>), py::arg("out"), py::arg("prec"));
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>
          (&refine_grid<D>), py::arg("out"), py::arg("inp"));
}
} // namespace vampyr
