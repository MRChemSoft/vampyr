#pragma once

#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/multiply.h>

namespace vampyr {
template <int D> void arithmetics(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("add",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            double,
                            FunctionTree<D> &,
                            double,
                            FunctionTree<D> &,
                            int,
                            bool>
                            (&mrcpp::add<D>),
          "prec"_a = -1.0,
          "out"_a,
          "a"_a = 1.0,
          "inp_a"_a,
          "b"_a = 1.0,
          "inp_b"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("add",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            std::vector<FunctionTree<D> *> &,
                            int,
                            bool>
                            (&mrcpp::add<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("add",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            std::vector<std::tuple<double, FunctionTree<D> *>>&,
                            int,
                            bool>
                            (&mrcpp::add<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("dot",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&mrcpp::dot<D>),
          "bra"_a,
          "ket"_a);

    m.def("multiply",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            double,
                            FunctionTree<D> &,
                            FunctionTree<D> &,
                            int,
                            bool,
                            bool>
                            (&mrcpp::multiply<D>),
          "prec"_a = -1.0,
          "out"_a,
          "c"_a = 1.0,
          "inp_a"_a,
          "inp_b"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def("multiply",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            std::vector<FunctionTree<D> *> &,
                            int,
                            bool,
                            bool>
                            (&mrcpp::multiply<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def("multiply",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            std::vector<std::tuple<double, FunctionTree<D> *>>&,
                            int,
                            bool,
                            bool>
                            (&mrcpp::multiply<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def("power",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            FunctionTree<D> &,
                            double,
                            int,
                            bool>
                            (&mrcpp::power<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "pow"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("square",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            FunctionTree<D> &,
                            int,
                            bool>
                            (&mrcpp::square<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
