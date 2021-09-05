#pragma once

#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/trees/TreeIterator.h>

namespace vampyr {
template <int D> void trees(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<MWTree<D>>(m, "MWTree")
        .def(py::init<MultiResolutionAnalysis<D>>(), "mra"_a)
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

    py::class_<FunctionTree<D>, MWTree<D>, RepresentableFunction<D>>(m, "FunctionTree")
        .def(py::init<MultiResolutionAnalysis<D>>())
        .def("getNGenNodes", &FunctionTree<D>::getNGenNodes)
        .def("deleteGenerated", &FunctionTree<D>::deleteGenerated)
        .def("integrate", &FunctionTree<D>::integrate)
        .def("normalize", &FunctionTree<D>::normalize)
        .def("rescale", &FunctionTree<D>::rescale, "coef"_a)
        .def("crop", &FunctionTree<D>::crop, "prec"_a, "split_fac"_a=1.0, "abs_prec"_a=false)
        .def("add", &FunctionTree<D>::add, "coef"_a=1.0, "inp"_a)
        .def("multiply", &FunctionTree<D>::multiply, "coef"_a=1.0, "inp"_a)
        .def("square", &FunctionTree<D>::square)
        .def("power", &FunctionTree<D>::power, "pow"_a)
        .def("saveTree", &FunctionTree<D>::saveTree, "filename"_a)
        .def("loadTree", &FunctionTree<D>::loadTree, "filename"_a);

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
        .def(py::init<Traverse, Iterator>(), "traverse"_a=TopDown, "iterator"_a=Lebesgue)
        .def(py::init<MWTree<D> &, Traverse, Iterator>(), "tree"_a, "traverse"_a=TopDown, "iterator"_a=Lebesgue)
        .def("setReturnGenNodes", &TreeIterator<D>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D>::setMaxDepth)
        .def("setTraverse", &TreeIterator<D>::setTraverse)
        .def("setIterator", &TreeIterator<D>::setIterator)
        .def("getNode", py::overload_cast<>(&TreeIterator<D>::getNode), py::return_value_policy::reference_internal)
        .def("init", &TreeIterator<D>::init)
        .def("next", &TreeIterator<D>::next);
}
} // namespace vampyr
