/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */


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

#include "operators/DerivativeOperator.h"
#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "treebuilders/add.h"
#include "treebuilders/multiply.h"
#include "treebuilders/apply.h"
#include "treebuilders/project.h"



using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void bases(py::module &);

PYBIND11_MODULE(vampyr2d, m) {

    m.doc() = "vampyr2d makes the MRCPP 2D functionality "
              "availible through a python interface";

    bases(m);

    m.def("project", py::overload_cast<double, FunctionTree<2>&,
        std::function<double (double, double)>, int>(&project2D),
        py::arg("precision"), py::arg("output_tree"), py::arg("function"), py::arg("maxIter")= -1);

    py::class_<PyBoundingBox<2>> (m, "BoundingBox")
        .def(py::init<int, py::array_t<int>, py::array_t <int>>(),
            "Initiates the BoundingBox of the area functions can exist")
        .def("getScale", &PyBoundingBox<2>::getScale);

    py::class_<MultiResolutionAnalysis<2>> (m, "MultiResolutionAnalysis")
        .def(py::init<PyBoundingBox<2>, ScalingBasis, int>(),
            py::arg("BoundingBox"), py::arg("Scaling basis"), py::arg("Scale"))
        .def("getOrder", &MultiResolutionAnalysis<2>::getOrder,
            "Returns the order of the scaling basis")
        .def("getMaxDepth", &MultiResolutionAnalysis<2>::getMaxDepth)
        .def("getMaxScale", &MultiResolutionAnalysis<2>::getMaxScale);

    py::class_<MWTree<2>> mwtree(m, "MWTree");
    mwtree
        .def(py::init<MultiResolutionAnalysis<2>>())
        .def("getSquareNorm", &MWTree<2>::getSquareNorm);

    py::class_<FunctionTree<2>> functree(m, "FunctionTree", mwtree);
    functree
        .def(py::init<MultiResolutionAnalysis<2>>())
        .def("integrate", &FunctionTree<2>::integrate,
            "Integral of the FunctionTree over all space")
        .def("clear", &FunctionTree<2>::clear,
            "Clears the FunctionTree")
        .def("normalize", &FunctionTree<2>::normalize,
            "Rescale the function by its norm")
        .def("rescale", &FunctionTree<2>::rescale,
            "Rescales the function")
        .def("evalf", py::overload_cast<double, double>(&FunctionTree<2>::evalf),
            "Returns the function value at a given point");

    m.def("apply", py::overload_cast<double, FunctionTree<2>&,
        ConvolutionOperator<2> &, FunctionTree<2> &, int>(&apply<2>),
        py::arg("precision"), py::arg("output_tree"), py::arg("oper"),
        py::arg("input_tree"), py::arg("maxIter") = -1);

    m.def("apply", py::overload_cast<FunctionTree<2>&, DerivativeOperator<2>&,
        FunctionTree<2> &, int>(&apply<2>),
        py::arg("out_tree"), py::arg("DerivativeOperator"),
        py::arg("input_tree"), py::arg("dir") = -1,
        "Applies a given derivative operator onto a FunctionTree");


    m.def("add", py::overload_cast<double, FunctionTree<2>&, double,
        FunctionTree<2> &, double, FunctionTree<2> &, int>(&add<2>),
        py::arg("precision"), py::arg("output_tree"), py::arg("a"),
        py::arg("tree_a"), py::arg("b"), py::arg("tree_b"),
        py::arg("maxIter") = -1,
        "Addition: output_tree = a*tree_a + b_tree_b");

    m.def("project", py::overload_cast<double, FunctionTree<2>&,
        RepresentableFunction<2> &, int>(&project<2>), py::arg("precision"),
        py::arg("output_tree"), py::arg("input_tree"), py::arg("maxIter")= -1);

    m.def("multiply", py::overload_cast<double, FunctionTree<2>&, double,
        FunctionTree<2> &, FunctionTree<2> &, int >(&multiply<2>),
        py::arg("precision"), py::arg("output_tree"), py::arg("c"), py::arg("tree_a"),
        py::arg("tree_b"), py::arg("maxIter") = -1,
        "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("dot", py::overload_cast<FunctionTree<2>&, FunctionTree<2>&>(&dot<2>),
        py::arg("bra"), py::arg("ket"));

    py::class_<DerivativeOperator<2>> deriv(m, "Derivative Operator");
    deriv
        .def(py::init<MultiResolutionAnalysis<2>>());

    py::class_<ABGVOperator<2>> (m, "ABGVOperator", deriv)
        .def(py::init< MultiResolutionAnalysis<2> &, double, double >(),
            "Derivative Operator: The ABGVOperator for differentiation");
}
} // namespace vampyr
