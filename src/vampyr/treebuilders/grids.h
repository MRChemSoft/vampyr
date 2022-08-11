#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/treebuilders/grid.h>

namespace vampyr {
template <int D> void advanced_grids(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("build_grid", py::overload_cast<FunctionTree<D> &, int>(&build_grid<D>), "out"_a, "scales"_a);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, const RepresentableFunction<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, std::vector<FunctionTree<D> *> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1);

    m.def(
        "build_grid",
        py::overload_cast<FunctionTree<D> &, std::vector<std::tuple<double, FunctionTree<D> *>> &, int>(&build_grid<D>),
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_grid<D>), "out"_a, "inp"_a);

    m.def("copy_func", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_func<D>), "out"_a, "inp"_a);

    m.def("clear_grid", py::overload_cast<FunctionTree<D> &>(&clear_grid<D>), "out"_a);

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, int>(&refine_grid<D>), "out"_a, "scales"_a);

    m.def("refine_grid",
          py::overload_cast<FunctionTree<D> &, double, bool>(&refine_grid<D>),
          "out"_a,
          "prec"_a,
          "abs_prec"_a = false);

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&refine_grid<D>), "out"_a, "inp"_a);
}

} // namespace vampyr
