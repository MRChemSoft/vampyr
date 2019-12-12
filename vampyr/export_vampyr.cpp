/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>

#include <string>

#include <MRCPP/trees/BandWidth.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>
#include <MRCPP/trees/OperatorTree.h>
#include <MRCPP/trees/SerialOperatorTree.h>
#include <MRCPP/trees/SerialTree.h>

#include "PyRepresentableFunction.h"
#include "applys/applys.h"
#include "functions/functions.h"
#include "math/math.h"
#include "project/project.h"
#include "trees/trees.h"
#include "world/world.h"

namespace py = pybind11;

namespace vampyr {
void bases(py::module &);
void operators(py::module &);
void operators_grid(py::module &);

template <int D> void bind_mr(py::module &mod) noexcept {
    applys<D>(mod);
    math<D>(mod);
    project<D>(mod);

    std::string name = "vampyr" + std::to_string(D) + "d";
    py::module sub_mod = mod.def_submodule(name.c_str());
    sub_mod.doc() = std::to_string(D) + "-dimensional bindings";

    operators(sub_mod);
    operators_grid(sub_mod);

    functions<D>(sub_mod);
    trees<D>(sub_mod);
    world<D>(sub_mod);
}

PYBIND11_MODULE(_vampyr, m) {
    using namespace mrcpp;
    m.doc() = "VAMPyR makes the MRCPP functionality available through a python interface";

    // Dimension-independent bindings go in the main module
    bases(m);

    py::class_<BandWidth> bandwidth(m, "BandWidth");
    bandwidth.def(py::init<int>())
        .def("clear", &BandWidth::clear)
        .def("isEmpty", &BandWidth::isEmpty)
        .def("getDepth", &BandWidth::getDepth)
        .def("getMaxWidth", &BandWidth::getMaxWidth)
        .def("setWidth", &BandWidth::setWidth);

    // NOTE What is the phantom_of_the_MWTree<2>? It's a dummy binding for
    // MWTree<2> to use exclusively to generate the bindings for its subclasses.
    py::class_<MWTree<2>> _mwtree(m, "phantom_of_the_MWTree2D");
    py::class_<OperatorTree>(m, "OperatorTree", _mwtree)
        .def(py::init<MultiResolutionAnalysis<2> &, double>())
        .def("getNormPrecision", &OperatorTree::getNormPrecision)
        .def("calcBandWidth", py::overload_cast<double>(&OperatorTree::calcBandWidth))
        .def("clearBandWidth", &OperatorTree::clearBandWidth)
        .def("setupOperNodeCache", &OperatorTree::setupOperNodeCache)
        .def("clearOperNodeCache", &OperatorTree::clearOperNodeCache)
        .def("mwTransformDown", &OperatorTree::mwTransformDown)
        .def("mwTransformUp", &OperatorTree::mwTransformUp);

    pybind11::class_<SerialTree<2>> _serialtree(m, "phantom_of_the_SerialTree2D");
    // <OperatorTree * wont work since raw pointers does not exist in python
    pybind11::class_<SerialOperatorTree>(m, "SerialOperatorTree", _serialtree)
        .def(py::init<OperatorTree *>())
        .def("allocRoots", &SerialOperatorTree::allocRoots)
        .def("allocChildren", &SerialOperatorTree::allocChildren)
        .def("allocGenChildren", &SerialOperatorTree::allocGenChildren)
        .def("deallocNodes", &SerialOperatorTree::deallocNodes)
        .def("deallocGenNodes", &SerialOperatorTree::deallocGenNodes)
        .def("deallocGenNodeChunks", &SerialOperatorTree::deallocGenNodeChunks);

    py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> _repfunc1d(
        m, "phantom_of_the_RepresentableFunction1D");
    legendre_poly(m, _repfunc1d);
    polynomial(m, _repfunc1d);
    boys_function(m, _repfunc1d);

    // Dimension-dependent bindings go into submodules
    //bind_mr<1>(m);
    //bind_mr<2>(m);
    bind_mr<3>(m);
}
} // namespace vampyr
