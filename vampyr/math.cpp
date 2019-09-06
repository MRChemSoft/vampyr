#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "trees/MWTree.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/OperatorTree.h"

#include "PyRepresentableFunction.h"

#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "operators/HelmholtzOperator.h"
#include "operators/PoissonOperator.h"
#include "operators/PHOperator.h"
#include "operators/BSOperator.h"
#include "operators/GreensKernel.h"
#include "operators/DerivativeKernel.h"
#include "operators/IdentityKernel.h"
#include "operators/MWOperator.h"
#include "operators/HelmholtzKernel.h"
#include "operators/OperatorState.h"

#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

#include "functions/GaussFunc.h"
#include "functions/GaussPoly.h"
#include "functions/Gaussian.h"
#include "functions/Polynomial.h"
#include "functions/RepresentableFunction.h"
//#include "functions/AnalyticFunction.h"
#include "functions/BoysFunction.h"
#include "functions/LegendrePoly.h"


using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void math(py::module &m) {
    const auto D = 3;

    m.def("add",
          py::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, double, FunctionTree<D> &, int>(
              &add<D>),
          "precision"_a,
          "output_tree"_a,
          "a"_a,
          "tree_a"_a,
          "b"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "Addition: output_tree = a*tree_a + b*tree_b");

    m.def("add",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int>(&add<D>),
          "precision"_a,
          "output_tree"_a,
          "FunctionTreeVector"_a,
          "maxIter"_a = -1,
          "Addition: output_tree = sum_i c_i f_i(x)");

    m.def("multiply",
          py::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, FunctionTree<D> &, int>(&multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "c"_a,
          "tree_a"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("multiply",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int>(&multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "inpit_vector"_a,
          "maxIter"_a = -1,
          "Multiplication: output_tree = prod_i c_i f_i(x)");

    m.def("dot", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&dot<D>), "bra"_a, "ket"_a);

    m.def("dot",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, FunctionTreeVector<D> &, int>(&dot<D>),
          "prec"_a,
          "out"_a,
          "inp_a"_a,
          "inp_b"_a,
          "maxIter"_a = -1);

    m.def("power",
          py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, double, int>(&power<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "power"_a,
          "maxIter"_a = -1);

    m.def("square",
          py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, int>(&square<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);
}
} // namespace vampyr
