#pragma once

#include <filesystem>

#include <pybind11/eigen.h>
#include <pybind11/stl/filesystem.h>

#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/trees/TreeIterator.h>

namespace vampyr {
template <int D>
auto impl__add__(mrcpp::FunctionTree<D> *inp_a, mrcpp::FunctionTree<D> *inp_b) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp_a->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({1.0, inp_b});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__sub__(mrcpp::FunctionTree<D> *inp_a, mrcpp::FunctionTree<D> *inp_b) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp_a->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({-1.0, inp_b});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__mul__(mrcpp::FunctionTree<D> *inp_a, mrcpp::FunctionTree<D> *inp_b) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp_a->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({1.0, inp_b});
    build_grid(*out, vec);
    build_grid(*out, 1);
    multiply(-1, *out, vec);
    return out;
};

template <int D>
auto impl__mul__(mrcpp::FunctionTree<D> *inp_a, double c) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp_a->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({c, inp_a});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__pos__(mrcpp::FunctionTree<D> *inp) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
    copy_grid(*out, *inp);
    copy_func(*out, *inp);
    return out;
};

template <int D>
auto impl__neg__(mrcpp::FunctionTree<D> *inp) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({-1.0, inp});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__truediv__(mrcpp::FunctionTree<D> *inp, double c) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
    FunctionTreeVector<D> vec;
    vec.push_back({1.0 / c, inp});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__pow__(mrcpp::FunctionTree<D> *inp, double c) -> std::unique_ptr<mrcpp::FunctionTree<D>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
    copy_grid(*out, *inp);
    copy_func(*out, *inp);
    refine_grid(*out, 1);
    out->power(c);
    return out;
};

template <int D> void trees(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<MWTree<D>>(m, "MWTree")
        .def("MRA", &MWTree<D>::getMRA, py::return_value_policy::reference_internal)
        .def("nNodes", &MWTree<D>::getNNodes)
        .def("nEndNodes", &MWTree<D>::getNEndNodes)
        .def("nRootNodes", &MWTree<D>::getNRootNodes)
        .def("fetchEndNode",
             py::overload_cast<int>(&MWTree<D>::getEndMWNode, py::const_),
             py::return_value_policy::reference_internal)
        .def("fetchRootNode",
             py::overload_cast<int>(&MWTree<D>::getRootMWNode, py::const_),
             py::return_value_policy::reference_internal)
        .def("rootScale", &MWTree<D>::getRootScale)
        .def("depth", &MWTree<D>::getDepth)
        .def("setZero",
             [](MWTree<D> *out) {
                 out->setZero();
                 return out;
             })
        .def("clear", &MWTree<D>::clear)
        .def("setName", &MWTree<D>::setName)
        .def("name", &MWTree<D>::getName)
        .def("fetchNode",
             py::overload_cast<NodeIndex<D>>(&MWTree<D>::getNode),
             py::return_value_policy::reference_internal)
        .def("squaredNorm", &MWTree<D>::getSquareNorm)
        .def("norm",
             [](MWTree<D> &tree) {
                 auto sqNorm = tree.getSquareNorm();
                 return (sqNorm >= 0.0) ? std::sqrt(sqNorm) : -1.0;
             })
        .def("__str__", [](MWTree<D> &tree) {
            std::ostringstream os;
            os << tree;
            return os.str();
        });

    py::class_<FunctionTree<D>, MWTree<D>, RepresentableFunction<D>>(m, "FunctionTree")
        .def(py::init<const MultiResolutionAnalysis<D> &, const std::string &>(), "mra"_a, "name"_a = "nn")
        .def("nGenNodes", &FunctionTree<D>::getNGenNodes)
        .def("deleteGenerated", &FunctionTree<D>::deleteGenerated)
        .def("integrate", &FunctionTree<D>::integrate)
        .def("normalize",
             [](FunctionTree<D> *out) {
                 out->normalize();
                 return out;
             })
        .def(
            "saveTree",
            [](FunctionTree<D> &obj, const std::string &filename) {
                namespace fs = std::filesystem;
                obj.saveTree(filename);
                return fs::absolute(fs::path(filename + ".tree"));
            },
            "filename"_a)
        .def("loadTree", &FunctionTree<D>::loadTree, "filename"_a)
        .def(
            "crop",
            [](FunctionTree<D> *out, double prec, bool abs_prec) {
                out->crop(prec, 1.0, abs_prec);
                return out;
            },
            "prec"_a,
            "abs_prec"_a = false)
        .def("deepCopy",
             [](FunctionTree<D> *inp) {
                 auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
                 copy_grid(*out, *inp);
                 copy_func(*out, *inp);
                 return out;
             })
        .def("quadrature",
             [](FunctionTree<D> *tree) {
                 if constexpr (D != 1) { throw std::runtime_error("quadrature only implemented for 1D"); }

                 // Current implementation only makes sense in 1D

                 std::vector<double> vec_pts;
                 // Iterate over all end nodes
                 for (int i = 0; i < tree->getNEndNodes(); i++) {
                     MWNode<D> &node = tree->getEndMWNode(i);

                     Eigen::MatrixXd pts;
                     node.getPrimitiveQuadPts(pts);

                     // Flatten the MatrixXd and add the points from this node to the vector
                     vec_pts.insert(vec_pts.end(), pts.data(), pts.data() + pts.size());
                 }

                 // Now we need to create an Eigen vector from our std::vector
                 Eigen::VectorXd final_pts =
                     Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(vec_pts.data(), vec_pts.size());

                 // Now final_pts holds all the points from all nodes
                 return final_pts;
             })
        .def("__call__", [](FunctionTree<D> &func, const Coord<D> &r) { return func.evalf_precise(r); })
        .def("__pos__", &impl__pos__<D>, py::is_operator())
        .def("__neg__", &impl__neg__<D>, py::is_operator())
        .def("__add__",  &impl__add__<D>, py::is_operator())
        .def("__iadd__", &impl__add__<D>, py::is_operator())
        .def("__sub__", &impl__sub__<D>, py::is_operator())
        .def("__isub__", &impl__sub__<D>, py::is_operator())
        .def("__mul__", py::overload_cast<FunctionTree<D> *, FunctionTree<D> *>(&impl__mul__<D>), py::is_operator())
        .def("__mul__", py::overload_cast<FunctionTree<D> *, double>(&impl__mul__<D>), py::is_operator())
        .def("__imul__", py::overload_cast<FunctionTree<D> *, FunctionTree<D> *>(&impl__mul__<D>), py::is_operator())
        .def("__imul__", py::overload_cast<FunctionTree<D> *, double>(&impl__mul__<D>), py::is_operator())
        .def("__rmul__", &impl__mul__<D>, py::is_operator())
        .def("__truediv__", &impl__truediv__<D>, py::is_operator())
        .def("__itruediv__", &impl__truediv__<D>, py::is_operator())
        .def("__pow__", &impl__pow__<D>, py::is_operator())
        .def("__ipow__", &impl__pow__<D>, py::is_operator());

    py::class_<MWNode<D>>(m, "MWNode")
        .def("depth", &MWNode<D>::getDepth)
        .def("scale", &MWNode<D>::getScale)
        .def("nCoefs", &MWNode<D>::getNCoefs)
        .def("nChildren", &MWNode<D>::getNChildren)
        .def("index",
             py::overload_cast<>(&MWNode<D>::getNodeIndex, py::const_),
             py::return_value_policy::reference_internal)
        .def("norm",
             [](MWNode<D> &node) {
                 auto sqNorm = node.getSquareNorm();
                 return (sqNorm >= 0.0) ? std::sqrt(sqNorm) : -1.0;
             })
        .def("squaredNorm", &MWNode<D>::getSquareNorm)
        .def("scalingNorm", &MWNode<D>::getScalingNorm)
        .def("waveletNorm", &MWNode<D>::getWaveletNorm)
        .def("componentNorm", &MWNode<D>::getComponentNorm)
        .def("isAllocated", &MWNode<D>::isAllocated)
        .def("isRootNode", &MWNode<D>::isRootNode)
        .def("isEndNode", &MWNode<D>::isEndNode)
        .def("isLeafNode", &MWNode<D>::isLeafNode)
        .def("isBranchNode", &MWNode<D>::isBranchNode)
        .def("isGenNode", &MWNode<D>::isGenNode)
        .def("hasParent", &MWNode<D>::hasParent)
        .def("hasCoefs", &MWNode<D>::hasCoefs)
        .def("quadrature",
             [](MWNode<D> &node) {
                 Eigen::MatrixXd pts;
                 node.getPrimitiveQuadPts(pts);
                 return pts;
             })
        .def("center", &MWNode<D>::getCenter)
        .def("upperBounds", &MWNode<D>::getUpperBounds)
        .def("lowerBounds", &MWNode<D>::getLowerBounds)
        .def("__str__", [](MWNode<D> &node) {
            std::ostringstream os;
            os << node;
            return os.str();
        });

    py::class_<FunctionNode<D>, MWNode<D>, std::unique_ptr<FunctionNode<D>, py::nodelete>>(m, "FunctionNode")
        .def("integrate", &FunctionNode<D>::integrate);

    py::class_<NodeIndex<D>>(m, "NodeIndex")
        .def(py::init<int, const std::array<int, D>>(), "scale"_a = 0, "translation"_a = std::array<int, D>{})
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("child", &NodeIndex<D>::child)
        .def("parent", &NodeIndex<D>::parent)
        .def("scale", &NodeIndex<D>::getScale)
        .def("setScale", &NodeIndex<D>::setScale)
        .def("translation", py::overload_cast<>(&NodeIndex<D>::getTranslation, py::const_))
        .def("translation", py::overload_cast<int>(&NodeIndex<D>::getTranslation, py::const_))
        .def("setTranslation", py::overload_cast<const std::array<int, D> &>(&NodeIndex<D>::setTranslation))
        .def("__str__", [](NodeIndex<D> &idx) {
            std::ostringstream os;
            os << idx;
            return os.str();
        });

    py::class_<TreeIterator<D>>(m, "TreeIterator")
        .def(py::init<Traverse, Iterator>(), "traverse"_a = TopDown, "iterator"_a = Lebesgue)
        .def(py::init<MWTree<D> &, Traverse, Iterator>(), "tree"_a, "traverse"_a = TopDown, "iterator"_a = Lebesgue)
        .def("setReturnGenNodes", &TreeIterator<D>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D>::setMaxDepth)
        .def("setTraverse", &TreeIterator<D>::setTraverse)
        .def("setIterator", &TreeIterator<D>::setIterator)
        .def("get", py::overload_cast<>(&TreeIterator<D>::getNode), py::return_value_policy::reference_internal)
        .def("init", &TreeIterator<D>::init)
        .def("next", &TreeIterator<D>::next);
}
} // namespace vampyr
