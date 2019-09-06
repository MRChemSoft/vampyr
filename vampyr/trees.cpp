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

void trees(py::module &m) {
    const auto D = 3;
    py::class_<MWTree<D>> mwtree(m, "MWTree");
    mwtree.def(py::init<MultiResolutionAnalysis<D>>()).def("getSquareNorm", &MWTree<D>::getSquareNorm);

    py::class_<FunctionTree<D>> functree(m, "FunctionTree", mwtree);
    functree.def(py::init<MultiResolutionAnalysis<D>>())
        .def("clear", &FunctionTree<D>::clear, "Clears the FunctionTree")
        .def("integrate", &FunctionTree<D>::integrate, "Integral of the FunctionTree over all space")
        .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")
        .def("rescale", &FunctionTree<D>::rescale, "Rescales the function")
        .def("getEndValues", &FunctionTree<D>::getEndValues, "")
        .def("saveTree", &FunctionTree<D>::saveTree, "Saves tree to file")
        .def("loadTree", &FunctionTree<D>::loadTree, "Loads tree to file")
        .def("crop", &FunctionTree<D>::crop, "Recurse down until an EndNode is found, and then crop children with too high precision")
        .def("multiply", &FunctionTree<D>::multiply, "Multiply the function by a given number and it is multiplied by the function")
        .def("add", &FunctionTree<D>::add, "Multiply the function by a given number and adds it to the function")
        .def("power", py::overload_cast<double>(&FunctionTree<D>::power), "Raise an existing function to a given power")
        .def("square", &FunctionTree<D>::square, "Multiply an existing function with itself")
        .def("evalf",
             py::overload_cast<const Coord<D> &>(&FunctionTree<D>::evalf, py::const_),
             "Returns the function value at a given point");


}
} // namespace vampyr
