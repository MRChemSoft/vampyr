#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/treebuilders/apply.h>

namespace vampyr {
template <int D> void applys(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // m.def("apply",
    //       py::overload_cast<double, FunctionTree<D> &, ConvolutionOperator<D> &, FunctionTree<D> &, int, bool>(
    //           &apply<D>),
    //       "precision"_a,
    //       "output_tree"_a,
    //       "oper"_a,
    //       "input_tree"_a,
    //       "maxIter"_a = -1,
    //       "abs_prec"_a = false,
    //       "Applies a given convolution operator onto a FunctionTree");

    m.def("apply",
          py::overload_cast<FunctionTree<D> &,
                            DerivativeOperator<D> &,
                            FunctionTree<D> &,
                            int>
                            (&apply<D>),
          "out"_a,
          "oper"_a,
          "inp"_a,
          "dir"_a = -1);

    m.def("divergence",
          py::overload_cast<FunctionTree<D> &,
                            DerivativeOperator<D> &,
                            std::vector<FunctionTree<D> *> &>
                            (&mrcpp::divergence<D>),
          "out"_a,
          "oper"_a,
          "inp"_a);

    m.def("gradient",
          py::overload_cast<DerivativeOperator<D> &,
                            FunctionTree<D> &>
                            (&gradient<D>),
          "oper"_a,
          "inp"_a);
}
} // namespace vampyr
