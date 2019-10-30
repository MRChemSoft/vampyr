#include <pybind11/pybind11.h>

#include "trees/FunctionTree.h"

#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"

#include "treebuilders/apply.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void applys(py::module &m) {
    m.def("apply",
          py::overload_cast<double, FunctionTree<3> &, ConvolutionOperator<3> &, FunctionTree<3> &, int>(&apply<3>),
          "precision"_a,
          "output_tree"_a,
          "oper"_a,
          "input_tree"_a,
          "maxIter"_a = -1,
          "Applies a given convolution operator onto a FunctionTree");

    m.def("apply",
          py::overload_cast<FunctionTree<3> &, DerivativeOperator<3> &, FunctionTree<3> &, int>(&apply<3>),
          "out_tree"_a,
          "DerivativeOperator"_a,
          "input_tree"_a,
          "dir"_a = -1,
          "Applies a given derivative operator onto a FunctionTree");

}
} // namespace vampyr
