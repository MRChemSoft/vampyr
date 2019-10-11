#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "trees/MultiResolutionAnalysis.h"

#include "PyRepresentableFunction.h"

#include "operators/DerivativeOperator.h"

#include "treebuilders/apply.h"
#include "treebuilders/grid.h"

#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void operators_grid(py::module &m) {

    const auto D = 3;
    
     m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, const RepresentableFunction<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTreeVector<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_grid<D>));
    m.def("copy_func", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_func<D>));

    m.def("clear_grid", py::overload_cast<FunctionTree<D> &>(&clear_grid<D>));
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, int>(&refine_grid<D>), "out"_a, "scales"_a);

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, double>(&refine_grid<D>), "out"_a, "prec"_a);
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&refine_grid<D>), "out"_a, "inp"_a);


}
} // namespace vampyr
