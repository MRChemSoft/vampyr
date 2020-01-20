#pragma once

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include <MRCPP/treebuilders/project.h>
#include <MRCPP/trees/FunctionTree.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D> void project(pybind11::module &m) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    m.def("project",
          pybind11::overload_cast<double, FunctionTree<D> &, std::function<double(const Coord<D> &)>, int, bool>(
              &mrcpp::project<D>),
          "precision"_a,
          "output_tree"_a,
          "function"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Projects an analytic function onto a FunctionTree");

    m.def("project",
          pybind11::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int, bool>(&mrcpp::project<D>),
          "precision"_a,
          "output_tree"_a,
          "GaussFunc"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
