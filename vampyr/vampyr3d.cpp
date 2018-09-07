/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <tuple>
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/eigen.h"
#include "pybind11/functional.h"

#include "trees/BoundingBox.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/MWTree.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "PyBoundingBox.h"
#include "PyAnalyticFunction.h"
#include "project.h"
#include "py_push_back.h"

#include "operators/DerivativeOperator.h"
#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/PoissonOperator.h"
#include "operators/HelmholtzOperator.h"

#include "treebuilders/add.h"
#include "treebuilders/multiply.h"
#include "treebuilders/apply.h"
#include "treebuilders/project.h"


using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void bases(py::module &);


PYBIND11_MODULE(vampyr3d, m) {

    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    bases(m);

    // FunctonTreeVector related start
    py::class_<FunctionTreeVector<3>> (m, "FunctionTreeVector")
        .def(py::init<>())
        .def("size", &FunctionTreeVector<3>::size);

    m.def("push_back", &py_push_back<3>);

    m.def("get_coef", &get_coef<3>);
    m.def("get_func", py::overload_cast<FunctionTreeVector<3> &, int>(&get_func<3>));
    // FunctonTreeVector related end
    m.def("project", py::overload_cast<double, FunctionTree<3>&,
        std::function<double (double, double, double)>, int>(&project3D),
        py::arg("precision"), py::arg("output_tree"), py::arg("function"), py::arg("maxIter")= -1,
        "Projects an analytic function onto a FunctionTree");

    py::class_<PyBoundingBox<3>> (m, "BoundingBox")
        .def(py::init<int, py::array_t<int>, py::array_t <int>>(),
            "Initiates the BoundingBox of the area functions can exist")
        .def("getScale", &PyBoundingBox<3>::getScale);

    py::class_<MultiResolutionAnalysis<3>> (m, "MultiResolutionAnalysis")
        .def(py::init<PyBoundingBox<3>, ScalingBasis, int>(),
            py::arg("BoundingBox"), py::arg("Basis Function"), py::arg("Scale"))
        .def("getOrder", &MultiResolutionAnalysis<3>::getOrder,
             "Returns the order of the scaling basis")
        .def("getMaxDepth", &MultiResolutionAnalysis<3>::getMaxDepth)
        .def("getMaxScale", &MultiResolutionAnalysis<3>::getMaxScale);

    py::class_<MWTree<3>> mwtree(m, "MWTree");
    mwtree
        .def(py::init<MultiResolutionAnalysis<3>>())
        .def("getSquareNorm", &MWTree<3>::getSquareNorm);

    py::class_<FunctionTree<3>> functree(m, "FunctionTree", mwtree);
    functree
        .def(py::init<MultiResolutionAnalysis<3>>())
        .def("integrate", &FunctionTree<3>::integrate,
            "Integral of the FunctionTree over all space")
        .def("clear", &FunctionTree<3>::clear,
            "Clears the FunctionTree")
        .def("normalize", &FunctionTree<3>::normalize,
            "Rescale the function by its norm, fixed grid")
        .def("rescale", &FunctionTree<3>::rescale,
            "Rescales the function")
        .def("evalf", py::overload_cast<double, double, double>(&FunctionTree<3>::evalf),
            "Returns the function value at a given point");


    m.def("apply", py::overload_cast<double, FunctionTree<3>&,
        ConvolutionOperator<3> &, FunctionTree<3> &, int>(&apply<3>),
        py::arg("precision"), py::arg("output_tree"), py::arg("oper"), py::arg("input_tree"),
        py::arg("maxIter") = -1,
        "Applies a given convolution operator onto a FunctionTree");

    m.def("apply", py::overload_cast<FunctionTree<3>&, DerivativeOperator<3>&,
        FunctionTree<3> &, int>(&apply<3>),
        py::arg("out_tree"), py::arg("DerivativeOperator"), py::arg("input_tree"),
        py::arg("dir") = -1,
        "Applies a given derivative operator onto a FunctionTree");

    m.def("add", py::overload_cast<double, FunctionTree<3>&, double,
        FunctionTree<3> &, double, FunctionTree<3> &, int>(&add<3>),
        py::arg("precision"), py::arg("output_tree"), py::arg("a"),  py::arg("tree_a"),
        py::arg("b"), py::arg("tree_b"), py::arg("maxIter") = -1,
        "Addition: output_tree = a*tree_a + b*tree_b");

    m.def("add", py::overload_cast<double, FunctionTree<3>&,
        FunctionTreeVector<3> &, int>(&add<3>),
        py::arg("precision"), py::arg("output_tree"), py::arg("FunctionTreeVector"), py::arg("maxIter") = -1,
        "Addition: output_tree = sum_i c_i f_i(x)");

    m.def("multiply", py::overload_cast<double, FunctionTree<3>&, double,
        FunctionTree<3> &, FunctionTree<3> &, int >(&multiply<3>),
        py::arg("precision"), py::arg("output_tree"), py::arg("c"), py::arg("tree_a"),
        py::arg("tree_b"), py::arg("maxIter") = -1,
        "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("multiply", py::overload_cast<double, FunctionTree<3>&, FunctionTreeVector<3> &, int >(&multiply<3>),
        py::arg("precision"), py::arg("output_tree"), py::arg("inpit_vector"),
        py::arg("maxIter") = -1,
        "Multiplication: output_tree = prod_i c_i f_i(x)");

    m.def("dot", py::overload_cast<FunctionTree<3>&, FunctionTree<3>&>(&dot<3>),
        py::arg("bra"), py::arg("ket"));

    py::class_<DerivativeOperator<3>> deriv(m, "Derivative Operator");
    deriv
        .def(py::init<MultiResolutionAnalysis<3>>());

    py::class_<ABGVOperator<3>> (m, "ABGVOperator", deriv)
        .def(py::init< MultiResolutionAnalysis<3> &, double, double>(),
            "Derivative Operator: The ABGVOperator for differentiation");

    py::class_<ConvolutionOperator<3>> convop(m, "ConvolutionOperator");
    convop
        .def(py::init<MultiResolutionAnalysis<3> &, double>());

    py::class_<PoissonOperator> (m, "PoissonOperator", convop)
        .def(py::init<MultiResolutionAnalysis<3> &, double >(),
             py::arg("MRA"), py::arg("precision"),
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");

    py::class_<HelmholtzOperator> (m, "HelmholtzOperator", convop)
        .def(py::init<MultiResolutionAnalysis<3> &, double, double>(),
            py::arg("MRA"), py::arg("mu"), py::arg("precision"),
            "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");
}
} // namespace vampyr
