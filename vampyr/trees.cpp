#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "trees/MWTree.h"
#include "trees/MultiResolutionAnalysis.h"

#include "PyRepresentableFunction.h"

#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "operators/HelmholtzOperator.h"
#include "operators/PoissonOperator.h"

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
        .def("integrate", &FunctionTree<D>::integrate, "Integral of the FunctionTree over all space")
        .def("clear", &FunctionTree<D>::clear, "Clears the FunctionTree")
        .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")
        .def("rescale", &FunctionTree<D>::rescale, "Rescales the function")
        .def("saveTree", &FunctionTree<D>::saveTree, "Saves tree to file")
        .def("loadTree", &FunctionTree<D>::loadTree, "Loads tree to file")
        .def("evalf",
             py::overload_cast<const Coord<D> &>(&FunctionTree<D>::evalf),
             "Returns the function value at a given point");

}
} // namespace vampyr
