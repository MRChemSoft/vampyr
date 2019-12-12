#include <pybind11/pybind11.h>

#include <MRCPP/operators/ConvolutionOperator.h>
#include <MRCPP/operators/DerivativeOperator.h>
#include <MRCPP/treebuilders/apply.h>
#include <MRCPP/trees/FunctionTree.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

template <int D> void applys(py::module &m) {
    m.def("apply",
          py::overload_cast<double,
                            mrcpp::FunctionTree<D> &,
                            mrcpp::ConvolutionOperator<D> &,
                            mrcpp::FunctionTree<D> &,
                            int,
                            bool>(&mrcpp::apply<D>),
          "precision"_a,
          "output_tree"_a,
          "oper"_a,
          "input_tree"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false,
          "Applies a given convolution operator onto a FunctionTree");

    m.def("apply",
          py::overload_cast<mrcpp::FunctionTree<D> &, mrcpp::DerivativeOperator<D> &, mrcpp::FunctionTree<D> &, int>(
              &mrcpp::apply<D>),
          "out_tree"_a,
          "DerivativeOperator"_a,
          "input_tree"_a,
          "dir"_a = -1,
          "Applies a given derivative operator onto a FunctionTree");
}
} // namespace vampyr
