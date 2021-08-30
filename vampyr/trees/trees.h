#pragma once

#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/trees/HilbertIterator.h>
#include <MRCPP/trees/LebesgueIterator.h>

namespace vampyr {
template <int D> void trees(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // The commented lines bellow were resulting in errors on compiling I (Evelin) didn't manage to fix.
    // I'm not completely sure of their importance on overall coding so left them to be "fixed" instead
    // of just deleting.
    py::class_<MWTree<D>>(m, "MWTree")
        .def(py::init<MultiResolutionAnalysis<D>>())
        .def("getMRA", &MWTree<D>::getMRA, py::return_value_policy::reference_internal)
        .def("getNNodes", &MWTree<D>::getNNodes)
        .def("getNEndNodes", &MWTree<D>::getNEndNodes)
        .def("getNRootNodes", &MWTree<D>::getNRootNodes)
        .def("getEndNode", py::overload_cast<int>(&MWTree<D>::getEndMWNode, py::const_), py::return_value_policy::reference_internal)
        .def("getRootNode", py::overload_cast<int>(&MWTree<D>::getRootMWNode, py::const_), py::return_value_policy::reference_internal)
        .def("getRootScale", &MWTree<D>::getRootScale)
        .def("getDepth", &MWTree<D>::getDepth)
        .def("setZero", &MWTree<D>::setZero)
        .def("clear", &MWTree<D>::clear)
        .def("setName", &MWTree<D>::setName)
        .def("getName", &MWTree<D>::getName)
        .def("getNode", py::overload_cast<NodeIndex<D>>(&MWTree<D>::getNode), py::return_value_policy::reference_internal)
        .def("getSquareNorm", &MWTree<D>::getSquareNorm)
        .def("__str__", [](MWTree<D> &tree) {
            std::ostringstream os;
            os << tree;
            return os.str();
        });
        // .def("makeNodeTable", py::overload_cast<MWNodeVector<D> &>(&MWTree<D>::makeNodeTable_))

    py::class_<FunctionTree<D>, MWTree<D>, RepresentableFunction<D>>(m, "FunctionTree")
        .def(py::init<MultiResolutionAnalysis<D>>())
        .def("getNGenNodes", &FunctionTree<D>::getNGenNodes)
        .def("integrate", &FunctionTree<D>::integrate)
        .def("deleteGenerated", &FunctionTree<D>::deleteGenerated)
        // .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")
        .def("saveTree", &FunctionTree<D>::saveTree, "filename"_a)
        .def("loadTree", &FunctionTree<D>::loadTree, "filename"_a);
        // .def("crop",
        //      &FunctionTree<D>::crop,
        //      "Recurse down until an EndNode is found, and then crop children with too high precision")
        // .def("multiply",
        //      &FunctionTree<D>::multiply,
        //      "Multiply the function by a given number and it is multiplied by the function")
        // .def("add", &FunctionTree<D>::add, "Multiply the function by a given number and adds it to the function")
        // .def("power", py::overload_cast<double>(&FunctionTree<D>::power), "Raise an existing function to a given power")
        // .def("square", &FunctionTree<D>::square, "Multiply an existing function with itself")
        // .def("rescale", py::overload_cast<double>(&FunctionTree<D>::rescale), "Rescales the function")

    py::class_<MWNode<D>>(m, "MWNode")
        .def("getDepth", &MWNode<D>::getDepth)
        .def("getScale", &MWNode<D>::getScale)
        .def("getNCoefs", &MWNode<D>::getNCoefs)
        .def("getNChildren", &MWNode<D>::getNChildren)
        .def("getNodeIndex", py::overload_cast<>(&MWNode<D>::getNodeIndex, py::const_), py::return_value_policy::reference_internal)
        .def("getSquareNorm", &MWNode<D>::getSquareNorm)
        .def("getScalingNorm", &MWNode<D>::getScalingNorm)
        .def("getWaveletNorm", &MWNode<D>::getWaveletNorm)
        .def("getComponentNorm", &MWNode<D>::getComponentNorm)
        .def("hasCoefs", &MWNode<D>::hasCoefs)
        .def("isAllocated", &MWNode<D>::isAllocated)
        .def("isRootNode", &MWNode<D>::isRootNode)
        .def("isEndNode", &MWNode<D>::isEndNode)
        .def("isLeafNode", &MWNode<D>::isLeafNode)
        .def("isBranchNode", &MWNode<D>::isBranchNode)
        .def("isGenNode", &MWNode<D>::isGenNode)
        .def("hasParent", &MWNode<D>::hasParent)
        .def("__str__", [](MWNode<D> &node) {
            std::ostringstream os;
            os << node;
            return os.str();
        });

    py::class_<FunctionNode<D>, MWNode<D>, std::unique_ptr<FunctionNode<D>, py::nodelete>>(m, "FunctionNode")
        .def("integrate", &FunctionNode<D>::integrate);

    py::class_<NodeIndex<D>>(m, "NodeIndex")
        .def(py::init<int, const std::array<int, D>>(),
             "scale"_a = 0,
             "translation"_a = std::array<int, D>{})
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("child", &NodeIndex<D>::child)
        .def("parent", &NodeIndex<D>::parent)
        .def("getScale", &NodeIndex<D>::getScale)
        .def("setScale", &NodeIndex<D>::setScale)
        .def("getTranslation", py::overload_cast<>(&NodeIndex<D>::getTranslation, py::const_))
        .def("getTranslation", py::overload_cast<int>(&NodeIndex<D>::getTranslation, py::const_))
        .def("setTranslation", py::overload_cast<const std::array<int, D> &>(&NodeIndex<D>::setTranslation))
        .def("__str__", [](NodeIndex<D> &idx) {
            std::ostringstream os;
            os << idx;
            return os.str();
        });

    py::class_<TreeIterator<D>>(m, "TreeIterator")
        .def("setReturnGenNodes", &TreeIterator<D>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D>::setMaxDepth)
        .def("getNode", py::overload_cast<>(&TreeIterator<D>::getNode), py::return_value_policy::reference_internal)
        .def("next", &TreeIterator<D>::next);

    py::class_<HilbertIterator<D>, TreeIterator<D>>(m, "HilbertIterator")
        .def(py::init<MWTree<D> *>(), "tree"_a);

    py::class_<LebesgueIterator<D>, TreeIterator<D>>(m, "LebesgueIterator")
        .def(py::init<MWTree<D> *>(), "tree"_a);
}
} // namespace vampyr
