#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/treebuilders/grid.h>

namespace vampyr {

template <int D> void grids(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, int scales) { build_grid<D, double>(out, scales); },
        "out"_a,
        "scales"_a);

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, FunctionTree<D, double> &inp, int max_iter) {
            build_grid<D, double>(out, inp, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);
}

template <int D> void advanced_grids(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, int scales) { build_grid<D, double>(out, scales); },
        "out"_a,
        "scales"_a);

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, FunctionTree<D, double> &inp, int max_iter) {
            build_grid<D, double>(out, inp, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, const RepresentableFunction<D, double> &inp, int max_iter) {
            build_grid<D, double>(out, inp, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, std::vector<FunctionTree<D, double> *> &inp, int max_iter) {
            FunctionTreeVector<D, double> vec;
            for (auto *tree : inp) vec.push_back({1.0, tree});
            build_grid<D, double>(out, vec, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "build_grid",
        [](FunctionTree<D, double> &out, std::vector<std::tuple<double, FunctionTree<D, double> *>> &inp, int max_iter) {
            FunctionTreeVector<D, double> vec;
            for (auto &t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
            build_grid<D, double>(out, vec, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "copy_grid",
        [](FunctionTree<D, double> &out, FunctionTree<D, double> &inp) { copy_grid<D, double>(out, inp); },
        "out"_a,
        "inp"_a);

    m.def(
        "copy_func",
        [](FunctionTree<D, double> &out, FunctionTree<D, double> &inp) { copy_func<D, double>(out, inp); },
        "out"_a,
        "inp"_a);

    m.def(
        "clear_grid", [](FunctionTree<D, double> &out) { clear_grid<D, double>(out); }, "out"_a);

    m.def(
        "refine_grid",
        [](FunctionTree<D, double> &out, int scales) { return refine_grid<D, double>(out, scales); },
        "out"_a,
        "scales"_a);

    m.def(
        "refine_grid",
        [](FunctionTree<D, double> &out, double prec, bool abs_prec) { return refine_grid<D, double>(out, prec, abs_prec); },
        "out"_a,
        "prec"_a,
        "abs_prec"_a = false);

    m.def(
        "refine_grid",
        [](FunctionTree<D, double> &out, FunctionTree<D, double> &inp) { return refine_grid<D, double>(out, inp); },
        "out"_a,
        "inp"_a);
}

} // namespace vampyr
