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

void project(py::module &m) {
    const auto D = 3;

    m.def(
        "project",
        py::overload_cast<double, FunctionTree<D> &, std::function<double(const Coord<D> &)>, int>(&mrcpp::project<D>),
        "precision"_a,
        "output_tree"_a,
        "function"_a,
        "maxIter"_a = -1,
        "Projects an analytic function onto a FunctionTree");

    m.def("project",
          py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int>(&mrcpp::project<D>),
          "precision"_a,
          "output_tree"_a,
          "GaussFunc"_a,
          "maxIter"_a = -1);
}
} // namespace vampyr
