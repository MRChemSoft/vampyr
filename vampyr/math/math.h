#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/multiply.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/FunctionTreeVector.h>

namespace vampyr {
template <int D> void math(pybind11::module &m) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    m.def("add",
          pybind11::
              overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, double, FunctionTree<D> &, int, bool>(
                  &add<D>),
          "precision"_a,
          "output_tree"_a,
          "a"_a,
          "tree_a"_a,
          "b"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Addition: output_tree = a*tree_a + b*tree_b");

    m.def("add",
          pybind11::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int, bool>(&add<D>),
          "precision"_a,
          "output_tree"_a,
          "FunctionTreeVector"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Addition: output_tree = sum_i c_i f_i(x)");

    m.def("multiply",
          pybind11::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, FunctionTree<D> &, int, bool>(
              &multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "c"_a,
          "tree_a"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("multiply",
          pybind11::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int, bool>(&multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "inpit_vector"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Multiplication: output_tree = prod_i c_i f_i(x)");

    m.def("dot", pybind11::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&dot<D>), "bra"_a, "ket"_a);

    m.def(
        "dot",
        pybind11::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, FunctionTreeVector<D> &, int, bool>(
            &dot<D>),
        "prec"_a,
        "out"_a,
        "inp_a"_a,
        "inp_b"_a,
        "maxIter"_a = -1,
        "abs_prec"_a = false);

    m.def("power",
          pybind11::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, double, int, bool>(&power<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "power"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false);

    m.def("square",
          pybind11::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, int, bool>(&square<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
