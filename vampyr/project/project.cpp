#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "trees/FunctionTree.h"

#include "PyRepresentableFunction.h"

#include "treebuilders/project.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void project(py::module &m) {
    const auto D = 3;

    m.def(
        "project",
        py::overload_cast<double, FunctionTree<D> &, std::function<double(const Coord<D> &)>, int, bool>(&mrcpp::project<D>),
        "precision"_a,
        "output_tree"_a,
        "function"_a,
        "maxIter"_a = -1,
        "abs_prec"_a = false,
        "Projects an analytic function onto a FunctionTree");

    m.def("project",
          py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int, bool>(&mrcpp::project<D>),
          "precision"_a,
          "output_tree"_a,
          "GaussFunc"_a,
          "maxIter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
