#pragma once

#include <pybind11/functional.h>

#include <MRCPP/treebuilders/project.h>

namespace vampyr {
template <int D> void project(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("project",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            std::function<double(const Coord<D> &r)>,
                            int,
                            bool>
                            (&mrcpp::project<D>),
          py::call_guard<py::gil_scoped_release>(),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("project",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            RepresentableFunction<D> &,
                            int,
                            bool>
                            (&mrcpp::project<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
