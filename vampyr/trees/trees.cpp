#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "trees/BoundingBox.h"
#include "trees/FunctionTree.h"
#include "trees/MWTree.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/OperatorTree.h"
#include "trees/NodeIndex.h"
#include "trees/MWNode.h"
#include "trees/BandWidth.h"
#include "trees/SerialTree.h"
#include "trees/SerialFunctionTree.h"
#include "trees/SerialOperatorTree.h"
#include "trees/TreeIterator.h"

#include "utils/mpi_utils.h"

#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void trees(py::module &m) {
    const auto D = 3;
   
    py::class_<BandWidth> bandwidth(m, "BandWidth");
    bandwidth.def(py::init<int>())
        .def("clear", &BandWidth::clear)
        .def("isEmpty", &BandWidth::isEmpty)
        .def("getDepth", &BandWidth::getDepth)
        .def("getMaxWidth", &BandWidth::getMaxWidth)
        .def("setWidth", &BandWidth::setWidth);


// The commented lines bellow were resulting in errors on compiling I (Evelin) didn't manage to fix.
// I'm not completely sure of their importance on overall coding so left them to be "fixed" instead 
// of just deleting.


    py::class_<MWTree<D>> mwtree(m, "MWTree");
    mwtree.def(py::init<MultiResolutionAnalysis<D>>())
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
        .def("mwTransform", py::overload_cast<int, bool>(&MWTree<D>::mwTransform))
        .def("setName", &MWTree<D>::setName)
        .def("getRootIndex", py::overload_cast<const Coord<D> &>(&MWTree<D>::getRootIndex, py::const_))
        //.def("makeNodeTable", py::overload_cast<MWNodeVector<D> &>(&MWTree<D>::makeNodeTable_))
        .def("resetEndNodeTable", &MWTree<D>::resetEndNodeTable)
        .def("clearEndNodeTable", &MWTree<D>::clearEndNodeTable)
        .def("deleteGenerated", &MWTree<D>::deleteGenerated)
        .def("getNThreads", &MWTree<D>::getNThreads)
        .def("saveTree", &MWTree<D>::saveTree)
        .def("loadTree", &MWTree<D>::loadTree)
        .def("countBranchNodes", py::overload_cast<int>(&MWTree<D>::countBranchNodes))
        .def("countLeafNodes", py::overload_cast<int>(&MWTree<D>::countLeafNodes))
        .def("countAllocNodes", py::overload_cast<int>(&MWTree<D>::countAllocNodes))
        .def("countNodes", py::overload_cast<int>(&MWTree<D>::countNodes))
        .def("RecountNodes", &MWTree<D>::RecountNodes);



    py::class_<FunctionTree<D>> functree(m, "FunctionTree", mwtree);
    functree.def(py::init<MultiResolutionAnalysis<D>>())
        .def("clear", &FunctionTree<D>::clear, "Clears the FunctionTree")
        .def("integrate", &FunctionTree<D>::integrate, "Integral of the FunctionTree over all space")
        .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")    
        .def("getEndValues", &FunctionTree<D>::getEndValues, "Returns the ending wavelength values")
        .def("setEndValues", &FunctionTree<D>::setEndValues, "Sets ending wavelength values")
        .def("saveTree", &FunctionTree<D>::saveTree, "Saves tree to file")
        .def("loadTree", &FunctionTree<D>::loadTree, "Loads tree to file")
        .def("crop", &FunctionTree<D>::crop, "Recurse down until an EndNode is found, and then crop children with too high precision")
        .def("multiply", &FunctionTree<D>::multiply, "Multiply the function by a given number and it is multiplied by the function")
        .def("add", &FunctionTree<D>::add, "Multiply the function by a given number and adds it to the function")
        .def("power", py::overload_cast<double>(&FunctionTree<D>::power), "Raise an existing function to a given power")
        .def("square", &FunctionTree<D>::square, "Multiply an existing function with itself")
        .def("rescale", py::overload_cast<double>(&FunctionTree<D>::rescale), "Rescales the function")
        .def("getNChunks", &FunctionTree<D>::getNChunks, "Get number of chunks")
        .def("getNChunksUsed", &FunctionTree<D>::getNChunksUsed, "Write size of tree")
        .def("printSerialIndices", &FunctionTree<D>::printSerialIndices, "Print Serial Indices")
        .def("evalf", py::overload_cast<const Coord<D> &>(&FunctionTree<D>::evalf, py::const_), "Returns the function value at a given point");



    py::class_<MWNode<D>> mwnode(m, "MWNode");
    mwnode.def(py::init<MWNode<D>>())
        .def("getSerialIx", &MWNode<D>::getSerialIx)
        .def("getTranslation", &MWNode<D>::getTranslation)
        .def("getCenter", &MWNode<D>::getCenter)
        .def("getBounds", &MWNode<D>::getBounds)
        .def("hasCoord", py::overload_cast<const Coord<D> &>(&MWNode<D>::hasCoord, py::const_))
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
       // .def("getCoefs", py::overload_cast<const double *>(&MWNode<D>::getCoefs, py::const_))
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


   // py::class_<GenNode<D>> gennode(m, "GenNode", funcnode);
    //gennode.def(py::init<FunctionNode<D>>());


    py::class_<HilbertPath<D>> hilbpath(m, "HilbertPath");
    hilbpath.def(py::init<HilbertPath<D>>())
        .def("getPath", &HilbertPath<D>::getPath)
        .def("getChildPath", &HilbertPath<D>::getChildPath)
        .def("getZIndex", &HilbertPath<D>::getZIndex)
        .def("getHIndex", &HilbertPath<D>::getHIndex);



    py::class_<NodeBox<D>> nodebox(m, "NodeBox");
    nodebox.def(py::init<BoundingBox<D>>())
       //.def("setNode", &NodeBox<D>::setNode);
       //.def("setNode", py::overload_cast<int, MWNode<D>**>(&NodeBox<D>::setNode))
        .def("clearNode", &NodeBox<D>::clearNode)      
        .def("getNode", py::overload_cast<const NodeIndex<D> &>(&NodeBox<D>::getNode, py::const_))
        .def("getNode", py::overload_cast<const Coord<D> &>(&NodeBox<D>::getNode, py::const_))
        .def("getNode", py::overload_cast<int>(&NodeBox<D>::getNode, py::const_))
        .def("getNOccupied", &NodeBox<D>::getNOccupied);


    // int * wont work since raw pointers does not exist in python
    py::class_<NodeIndex<D>> nodeindex(m, "NodeIndex");
    nodeindex.def(py::init<int, int *>())  
        .def("setScale", &NodeIndex<D>::setScale)
        .def("setTranslation", &NodeIndex<D>::setTranslation)
        .def("getScale", &NodeIndex<D>::getScale)
        .def("getTranslation", py::overload_cast<int>(&NodeIndex<D>::getTranslation, py::const_));
      //  .def("getTranslation", py::overload_cast<int *>(&NodeIndex<D>::getTranslation))
      //  .def("getTranslation", py::overload_cast<const int *>(&NodeIndex<D>::getTranslation, py::const_));



    py::class_<OperatorTree>(m, "OperatorTree", mwtree)
        .def(py::init<MultiResolutionAnalysis<2> &, double>())      
        .def("getNormPrecision", &OperatorTree::getNormPrecision)
        .def("calcBandWidth", py::overload_cast<double>(&OperatorTree::calcBandWidth))
        .def("clearBandWidth", &OperatorTree::clearBandWidth)
        .def("setupOperNodeCache", &OperatorTree::setupOperNodeCache)
        .def("clearOperNodeCache", &OperatorTree::clearOperNodeCache)
        .def("mwTransformDown", &OperatorTree::mwTransformDown)
        .def("mwTransformUp", &OperatorTree::mwTransformUp);



    py::class_<SerialTree<D>> serialtree(m, "SerialTree");
    serialtree.def("isShared", &SerialTree<D>::isShared)
        .def("S_mwTransform", &SerialTree<D>::S_mwTransform)
        .def("S_mwTransformBack", &SerialTree<D>::S_mwTransformBack);


    // <FunctionTree<D> * and SharedMemory * wont work since raw pointers does not exist in python
    py::class_<SerialFunctionTree<D>> serialfunctree(m, "SerialFunctionTree", serialtree);
    serialfunctree.def(py::init<FunctionTree<D> *, SharedMemory *>())   
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


    // <OperatorTree * wont work since raw pointers does not exist in python
    py::class_<SerialOperatorTree>(m, "SerialOperatorTree", serialtree)
        .def(py::init<OperatorTree *>())  
        .def("allocRoots", &SerialOperatorTree::allocRoots)
        .def("allocChildren", &SerialOperatorTree::allocChildren)
        .def("allocGenChildren", &SerialOperatorTree::allocGenChildren)
        .def("deallocNodes", &SerialOperatorTree::deallocNodes)
        .def("deallocGenNodes", &SerialOperatorTree::deallocGenNodes)
        .def("deallocGenNodeChunks", &SerialOperatorTree::deallocGenNodeChunks);


    py::class_<TreeIterator<D>> treeiterator(m, "TreeIterator");
    treeiterator.def("setReturnGenNodes", &TreeIterator<D>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D>::setMaxDepth)
        .def("next", &TreeIterator<D>::next);


}
} // namespace vampyr
