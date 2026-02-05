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
auto impl__add__(mrcpp::FunctionTree<D, double> *inp_a, mrcpp::FunctionTree<D, double> *inp_b)
    -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp_a->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({1.0, inp_b});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__sub__(mrcpp::FunctionTree<D, double> *inp_a, mrcpp::FunctionTree<D, double> *inp_b)
    -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp_a->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({-1.0, inp_b});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__mul__(mrcpp::FunctionTree<D, double> *inp_a, mrcpp::FunctionTree<D, double> *inp_b)
    -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp_a->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({1.0, inp_a});
    vec.push_back({1.0, inp_b});
    build_grid(*out, vec);
    build_grid(*out, 1);
    multiply(-1, *out, vec);
    return out;
};

template <int D>
auto impl__mul__(mrcpp::FunctionTree<D, double> *inp_a, double c) -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp_a->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({c, inp_a});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D> auto impl__pos__(mrcpp::FunctionTree<D, double> *inp) -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
    copy_grid(*out, *inp);
    copy_func(*out, *inp);
    return out;
};

template <int D> auto impl__neg__(mrcpp::FunctionTree<D, double> *inp) -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({-1.0, inp});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D>
auto impl__truediv__(mrcpp::FunctionTree<D, double> *inp, double c) -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
    FunctionTreeVector<D, double> vec;
    vec.push_back({1.0 / c, inp});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
};

template <int D> auto impl__pow__(mrcpp::FunctionTree<D, double> *inp, double c) -> std::unique_ptr<mrcpp::FunctionTree<D, double>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
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

    py::class_<MWTree<D, double>>(m, "MWTree")
        .def("MRA", &MWTree<D, double>::getMRA, py::return_value_policy::reference_internal)
        .def("nNodes", &MWTree<D, double>::getNNodes)
        .def("nEndNodes", &MWTree<D, double>::getNEndNodes)
        .def("nRootNodes", &MWTree<D, double>::getNRootNodes)
        .def("fetchEndNode",
             [](MWTree<D, double>& tree, int i) -> MWNode<D, double>& { return tree.getEndMWNode(i); },
             py::return_value_policy::reference_internal)
        .def("fetchRootNode",
             [](MWTree<D, double>& tree, int i) -> MWNode<D, double>& { return tree.getRootMWNode(i); },
             py::return_value_policy::reference_internal)
        .def("rootScale", &MWTree<D, double>::getRootScale)
        .def("depth", &MWTree<D, double>::getDepth)
        .def("setZero",
             [](MWTree<D, double> *out) {
                 out->setZero();
                 return out;
             })
        .def("clear", &MWTree<D, double>::clear)
        .def("setName", &MWTree<D, double>::setName)
        .def("name", &MWTree<D, double>::getName)
        .def("fetchNode",
             [](MWTree<D, double>& tree, NodeIndex<D> idx) -> MWNode<D, double>& { return tree.getNode(idx); },
             py::return_value_policy::reference_internal)
        .def("squaredNorm", &MWTree<D, double>::getSquareNorm)
        .def("norm",
             [](MWTree<D, double> &tree) {
                 auto sqNorm = tree.getSquareNorm();
                 return (sqNorm >= 0.0) ? std::sqrt(sqNorm) : -1.0;
             })
        .def("__str__", [](MWTree<D, double> &tree) {
            std::ostringstream os;
            os << tree;
            return os.str();
        });

    py::class_<FunctionTree<D, double>, MWTree<D, double>, RepresentableFunction<D, double>>(m, "FunctionTree")
        .def(py::init<const MultiResolutionAnalysis<D> &, const std::string &>(), "mra"_a, "name"_a = "nn")
        .def("nGenNodes", &FunctionTree<D, double>::getNGenNodes)
        .def("deleteGenerated", &FunctionTree<D, double>::deleteGenerated)
        .def("integrate", &FunctionTree<D, double>::integrate)
        .def("normalize",
             [](FunctionTree<D, double> *out) {
                 out->normalize();
                 return out;
             })
        .def(
            "saveTree",
            [](FunctionTree<D, double> &obj, const std::string &filename) {
                namespace fs = std::filesystem;
                obj.saveTree(filename);
                return fs::absolute(fs::path(filename + ".tree"));
            },
            "filename"_a)
        .def("loadTree", &FunctionTree<D, double>::loadTree, "filename"_a)
        .def(
            "crop",
            [](FunctionTree<D, double> *out, double prec, bool abs_prec) {
                out->crop(prec, 1.0, abs_prec);
                return out;
            },
            "prec"_a,
            "abs_prec"_a = false)
        .def("deepCopy",
             [](FunctionTree<D, double> *inp) {
                 auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
                 copy_grid(*out, *inp);
                 copy_func(*out, *inp);
                 return out;
             })
        .def("__call__", [](FunctionTree<D, double> &func, const Coord<D> &r) { return func.evalf_precise(r); })
        .def("__pos__", &impl__pos__<D>, py::is_operator())
        .def("__neg__", &impl__neg__<D>, py::is_operator())
        .def("__add__", &impl__add__<D>, py::is_operator())
        .def("__iadd__", &impl__add__<D>, py::is_operator())
        .def("__sub__", &impl__sub__<D>, py::is_operator())
        .def("__isub__", &impl__sub__<D>, py::is_operator())
        .def("__mul__", py::overload_cast<FunctionTree<D, double> *, FunctionTree<D, double> *>(&impl__mul__<D>), py::is_operator())
        .def("__mul__", py::overload_cast<FunctionTree<D, double> *, double>(&impl__mul__<D>), py::is_operator())
        .def("__imul__", py::overload_cast<FunctionTree<D, double> *, FunctionTree<D, double> *>(&impl__mul__<D>), py::is_operator())
        .def("__imul__", py::overload_cast<FunctionTree<D, double> *, double>(&impl__mul__<D>), py::is_operator())
        .def("__rmul__", py::overload_cast<FunctionTree<D, double> *, double>(&impl__mul__<D>), py::is_operator())
        .def("__truediv__", &impl__truediv__<D>, py::is_operator())
        .def("__itruediv__", &impl__truediv__<D>, py::is_operator())
        .def("__pow__", &impl__pow__<D>, py::is_operator())
        .def("__ipow__", &impl__pow__<D>, py::is_operator());

    py::class_<MWNode<D, double>>(m, "MWNode")
        .def("depth", &MWNode<D, double>::getDepth)
        .def("scale", &MWNode<D, double>::getScale)
        .def("nCoefs", &MWNode<D, double>::getNCoefs)
        .def("nChildren", &MWNode<D, double>::getNChildren)
        .def("index",
             [](MWNode<D, double>& node) -> const NodeIndex<D>& { return node.getNodeIndex(); },
             py::return_value_policy::reference_internal)
        .def("norm",
             [](MWNode<D, double> &node) {
                 auto sqNorm = node.getSquareNorm();
                 return (sqNorm >= 0.0) ? std::sqrt(sqNorm) : -1.0;
             })
        .def("squaredNorm", &MWNode<D, double>::getSquareNorm)
        .def("scalingNorm", &MWNode<D, double>::getScalingNorm)
        .def("waveletNorm", &MWNode<D, double>::getWaveletNorm)
        .def("componentNorm", &MWNode<D, double>::getComponentNorm)
        .def("isAllocated", &MWNode<D, double>::isAllocated)
        .def("isRootNode", &MWNode<D, double>::isRootNode)
        .def("isEndNode", &MWNode<D, double>::isEndNode)
        .def("isLeafNode", &MWNode<D, double>::isLeafNode)
        .def("isBranchNode", &MWNode<D, double>::isBranchNode)
        .def("isGenNode", &MWNode<D, double>::isGenNode)
        .def("hasParent", &MWNode<D, double>::hasParent)
        .def("hasCoefs", &MWNode<D, double>::hasCoefs)
        .def("center", &MWNode<D, double>::getCenter)
        .def("upperBounds", &MWNode<D, double>::getUpperBounds)
        .def("lowerBounds", &MWNode<D, double>::getLowerBounds)
        .def("__str__", [](MWNode<D, double> &node) {
            std::ostringstream os;
            os << node;
            return os.str();
        });

    py::class_<FunctionNode<D, double>, MWNode<D, double>, std::unique_ptr<FunctionNode<D, double>, py::nodelete>>(m, "FunctionNode")
        .def("integrate", &FunctionNode<D, double>::integrate);

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

    py::class_<TreeIterator<D, double>>(m, "TreeIterator")
        .def(py::init<Traverse, Iterator>(), "traverse"_a = TopDown, "iterator"_a = Lebesgue)
        .def(py::init<MWTree<D, double> &, Traverse, Iterator>(), "tree"_a, "traverse"_a = TopDown, "iterator"_a = Lebesgue)
        .def("setReturnGenNodes", &TreeIterator<D, double>::setReturnGenNodes)
        .def("setMaxDepth", &TreeIterator<D, double>::setMaxDepth)
        .def("setTraverse", &TreeIterator<D, double>::setTraverse)
        .def("setIterator", &TreeIterator<D, double>::setIterator)
        .def("get", [](TreeIterator<D, double>& it) -> MWNode<D, double>& { return it.getNode(); }, py::return_value_policy::reference_internal)
        .def("init", &TreeIterator<D, double>::init)
        .def("next", &TreeIterator<D, double>::next);
}
} // namespace vampyr
