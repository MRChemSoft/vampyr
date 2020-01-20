#pragma once

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/apply.h>
#include <MRCPP/treebuilders/grid.h>
#include <MRCPP/treebuilders/multiply.h>
#include <MRCPP/treebuilders/project.h>
#include <MRCPP/trees/BoundingBox.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>
#include <MRCPP/trees/NodeIndex.h>
#include <MRCPP/trees/SerialFunctionTree.h>
#include <MRCPP/trees/SerialTree.h>
#include <MRCPP/trees/TreeIterator.h>
#include <MRCPP/utils/mpi_utils.h>

namespace vampyr {
template <int D> void trees(pybind11::module &m) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    // The commented lines bellow were resulting in errors on compiling I (Evelin) didn't manage to fix.
    // I'm not completely sure of their importance on overall coding so left them to be "fixed" instead
    // of just deleting.
    pybind11::class_<MWTree<D>> mwtree(m, "MWTree");
    mwtree.def(pybind11::init<MultiResolutionAnalysis<D>>())
        .def("getSquareNorm", &MWTree<D>::getSquareNorm)
        .def("calcSquareNorm", &MWTree<D>::calcSquareNorm)
        .def("clearSquareNorm", &MWTree<D>::clearSquareNorm)
        .def("getOrder", &MWTree<D>::getOrder)
        .def("getKp1", &MWTree<D>::getKp1)
        .def("getKp1_d", &MWTree<D>::getKp1_d)
        .def("getDim", &MWTree<D>::getDim)
        .def("getTDim", &MWTree<D>::getTDim)
        .def("getNNodes", &MWTree<D>::getNNodes)
        .def("getNEndNodes", &MWTree<D>::getNEndNodes)
        .def("getNGenNodes", &MWTree<D>::getNGenNodes)
        .def("getRootScale", &MWTree<D>::getRootScale)
        .def("getDepth", &MWTree<D>::getDepth)
        .def("getSizeNodes", &MWTree<D>::getSizeNodes)
        .def("mwTransform", pybind11::overload_cast<int, bool>(&MWTree<D>::mwTransform))
        .def("setName", &MWTree<D>::setName)
        .def("getRootIndex", pybind11::overload_cast<const Coord<D> &>(&MWTree<D>::getRootIndex, pybind11::const_))
        //.def("makeNodeTable", pybind11::overload_cast<MWNodeVector<D> &>(&MWTree<D>::makeNodeTable_))
        .def("resetEndNodeTable", &MWTree<D>::resetEndNodeTable)
        .def("clearEndNodeTable", &MWTree<D>::clearEndNodeTable)
        .def("deleteGenerated", &MWTree<D>::deleteGenerated)
        .def("getNThreads", &MWTree<D>::getNThreads)
        .def("saveTree", &MWTree<D>::saveTree)
        .def("loadTree", &MWTree<D>::loadTree)
        .def("countBranchNodes", pybind11::overload_cast<int>(&MWTree<D>::countBranchNodes))
        .def("countLeafNodes", pybind11::overload_cast<int>(&MWTree<D>::countLeafNodes))
        .def("countAllocNodes", pybind11::overload_cast<int>(&MWTree<D>::countAllocNodes))
        .def("countNodes", pybind11::overload_cast<int>(&MWTree<D>::countNodes))
        .def("RecountNodes", &MWTree<D>::RecountNodes);

    pybind11::class_<FunctionTree<D>> functree(m, "FunctionTree", mwtree);
    functree.def(pybind11::init<MultiResolutionAnalysis<D>>())
        .def("clear", &FunctionTree<D>::clear, "Clears the FunctionTree")
        .def("integrate", &FunctionTree<D>::integrate, "Integral of the FunctionTree over all space")
        .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")
        .def("getEndValues", &FunctionTree<D>::getEndValues, "Returns the ending wavelength values")
        .def("setEndValues", &FunctionTree<D>::setEndValues, "Sets ending wavelength values")
        .def("saveTree", &FunctionTree<D>::saveTree, "Saves tree to file")
        .def("loadTree", &FunctionTree<D>::loadTree, "Loads tree to file")
        .def("crop",
             &FunctionTree<D>::crop,
             "Recurse down until an EndNode is found, and then crop children with too high precision")
        .def("multiply",
             &FunctionTree<D>::multiply,
             "Multiply the function by a given number and it is multiplied by the function")
        .def("add", &FunctionTree<D>::add, "Multiply the function by a given number and adds it to the function")
        .def("power", pybind11::overload_cast<double>(&FunctionTree<D>::power), "Raise an existing function to a given power")
        .def("square", &FunctionTree<D>::square, "Multiply an existing function with itself")
        .def("rescale", pybind11::overload_cast<double>(&FunctionTree<D>::rescale), "Rescales the function")
        .def("getNChunks", &FunctionTree<D>::getNChunks, "Get number of chunks")
        .def("getNChunksUsed", &FunctionTree<D>::getNChunksUsed, "Write size of tree")
        .def("printSerialIndices", &FunctionTree<D>::printSerialIndices, "Print Serial Indices")
        .def("evalf",
             pybind11::overload_cast<const Coord<D> &>(&FunctionTree<D>::evalf, pybind11::const_),
             "Returns the function value at a given point");

    pybind11::class_<MWNode<D>> mwnode(m, "MWNode");
    mwnode.def(pybind11::init<MWNode<D>>())
        .def("getSerialIx", &MWNode<D>::getSerialIx)
        .def("getTranslation", &MWNode<D>::getTranslation)
        .def("getCenter", &MWNode<D>::getCenter)
        .def("getBounds", &MWNode<D>::getBounds)
        .def("hasCoord", pybind11::overload_cast<const Coord<D> &>(&MWNode<D>::hasCoord, pybind11::const_))
        .def("isCompatible", &MWNode<D>::isCompatible)
        .def("isAncestor", &MWNode<D>::isAncestor)
        .def("isDecendant", &MWNode<D>::isDecendant)
        .def("hasCoefs", &MWNode<D>::hasCoefs)
        .def("isRootNode", &MWNode<D>::isRootNode)
        .def("isEndNode", &MWNode<D>::isEndNode)
        .def("isGenNode", &MWNode<D>::isGenNode)
        .def("isLeafNode", &MWNode<D>::isLeafNode)
        .def("isAllocated", &MWNode<D>::isAllocated)
        .def("isBranchNode", &MWNode<D>::isBranchNode)
        .def("isLooseNode", &MWNode<D>::isLooseNode)
        .def("getSquareNorm", &MWNode<D>::getSquareNorm)
        .def("getScalingNorm", &MWNode<D>::getScalingNorm)
        .def("getWaveletNorm", &MWNode<D>::getWaveletNorm)
        .def("getComponentNorm", &MWNode<D>::getComponentNorm)
        //.def("hasComponentNorms", &MWNode<D>::hasComponentNorms)
        .def("getNCoefs", &MWNode<D>::getNCoefs)
        // .def("getCoefs", pybind11::overload_cast<const double *>(&MWNode<D>::getCoefs, pybind11::const_))
        .def("printCoefs", &MWNode<D>::printCoefs)
        .def("zeroCoefs", &MWNode<D>::zeroCoefs)
        .def("setCoefBlock", &MWNode<D>::setCoefBlock)
        .def("addCoefBlock", &MWNode<D>::addCoefBlock)
        .def("zeroCoefBlock", &MWNode<D>::zeroCoefBlock)
        .def("calcNorms", &MWNode<D>::calcNorms)
        .def("zeroNorms", &MWNode<D>::zeroNorms)
        .def("clearNorms", &MWNode<D>::clearNorms)
        .def("createChildren", &MWNode<D>::createChildren)
        .def("genChildren", &MWNode<D>::genChildren)
        .def("deleteChildren", &MWNode<D>::deleteChildren)
        .def("cvTransform", &MWNode<D>::cvTransform)
        .def("mwTransform", &MWNode<D>::mwTransform)
        .def("splitCheck", &MWNode<D>::splitCheck)
        .def("setHasCoefs", &MWNode<D>::setHasCoefs)
        .def("setIsEndNode", &MWNode<D>::setIsEndNode)
        .def("setIsGenNode", &MWNode<D>::setIsGenNode)
        .def("setIsRootNode", &MWNode<D>::setIsRootNode)
        .def("setIsLeafNode", &MWNode<D>::setIsLeafNode)
        .def("setIsAllocated", &MWNode<D>::setIsAllocated)
        .def("setIsBranchNode", &MWNode<D>::setIsBranchNode)
        .def("setIsLooseNode", &MWNode<D>::setIsLooseNode)
        .def("clearHasCoefs", &MWNode<D>::clearHasCoefs)
        .def("clearIsEndNode", &MWNode<D>::clearIsEndNode)
        .def("clearIsGenNode", &MWNode<D>::clearIsGenNode)
        .def("clearIsRootNode", &MWNode<D>::clearIsRootNode)
        .def("clearIsAllocated", &MWNode<D>::clearIsAllocated);

    // pybind11::class_<GenNode<D>> gennode(m, "GenNode", funcnode);
    //gennode.def(pybind11::init<FunctionNode<D>>());

    pybind11::class_<HilbertPath<D>> hilbpath(m, "HilbertPath");
    hilbpath.def(pybind11::init<HilbertPath<D>>())
        .def("getPath", &HilbertPath<D>::getPath)
        .def("getChildPath", &HilbertPath<D>::getChildPath)
        .def("getZIndex", &HilbertPath<D>::getZIndex)
        .def("getHIndex", &HilbertPath<D>::getHIndex);

    pybind11::class_<NodeBox<D>> nodebox(m, "NodeBox");
    nodebox
        .def(pybind11::init<BoundingBox<D>>())
        //.def("setNode", &NodeBox<D>::setNode);
        //.def("setNode", pybind11::overload_cast<int, MWNode<D>**>(&NodeBox<D>::setNode))
        .def("clearNode", &NodeBox<D>::clearNode)
        .def("getNode", pybind11::overload_cast<const NodeIndex<D> &>(&NodeBox<D>::getNode, pybind11::const_))
        .def("getNode", pybind11::overload_cast<const Coord<D> &>(&NodeBox<D>::getNode, pybind11::const_))
        .def("getNode", pybind11::overload_cast<int>(&NodeBox<D>::getNode, pybind11::const_))
        .def("getNOccupied", &NodeBox<D>::getNOccupied);

    // int * wont work since raw pointers does not exist in python
    pybind11::class_<NodeIndex<D>> nodeindex(m, "NodeIndex");
    nodeindex.def(pybind11::init<int, int *>())
        .def("setScale", &NodeIndex<D>::setScale)
        .def("setTranslation", &NodeIndex<D>::setTranslation)
        .def("getScale", &NodeIndex<D>::getScale)
        .def("getTranslation", pybind11::overload_cast<int>(&NodeIndex<D>::getTranslation, pybind11::const_));
    //  .def("getTranslation", pybind11::overload_cast<int *>(&NodeIndex<D>::getTranslation))
    //  .def("getTranslation", pybind11::overload_cast<const int *>(&NodeIndex<D>::getTranslation, pybind11::const_));

    pybind11::class_<SerialTree<D>> serialtree(m, "SerialTree");
    serialtree.def("isShared", &SerialTree<D>::isShared)
        .def("S_mwTransform", &SerialTree<D>::S_mwTransform)
        .def("S_mwTransformBack", &SerialTree<D>::S_mwTransformBack);

    // <FunctionTree<D> * and SharedMemory * wont work since raw pointers does not exist in python
    pybind11::class_<SerialFunctionTree<D>> serialfunctree(m, "SerialFunctionTree", serialtree);
    serialfunctree.def(pybind11::init<FunctionTree<D> *, SharedMemory *>())
        .def("allocRoots", &SerialFunctionTree<D>::allocRoots)
        .def("allocChildren", &SerialFunctionTree<D>::allocChildren)
        .def("allocGenChildren", &SerialFunctionTree<D>::allocGenChildren)
        .def("deallocNodes", &SerialFunctionTree<D>::deallocNodes)
        .def("deallocGenNodes", &SerialFunctionTree<D>::deallocGenNodes)
        .def("deallocGenNodeChunks", &SerialFunctionTree<D>::deallocGenNodeChunks)
        .def("getNChunks", &SerialFunctionTree<D>::getNChunks)
        .def("getNChunksUsed", &SerialFunctionTree<D>::getNChunksUsed)
        .def("shrinkChunks", &SerialFunctionTree<D>::shrinkChunks)
        .def("rewritePointers", &SerialFunctionTree<D>::rewritePointers)
        .def("clear", &SerialFunctionTree<D>::clear);

    pybind11::class_<TreeIterator<D>> treeiterator(m, "TreeIterator");
    treeiterator.def("setReturnGenNodes", &TreeIterator<D>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D>::setMaxDepth)
        .def("next", &TreeIterator<D>::next);
}
} // namespace vampyr
