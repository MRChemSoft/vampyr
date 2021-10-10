#pragma once

#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/multiply.h>

namespace vampyr {
template <int D> void arithmetics(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("sum", [] (std::vector<FunctionTree<D> *> &inp) {
              auto out = std::unique_ptr<FunctionTree<D>>(nullptr);
              if (inp.size() > 0) {
                  auto &mra = inp[0]->getMRA();
                  out = std::make_unique<FunctionTree<D>>(mra);
                  build_grid(*out, inp);
                  add(-1.0, *out, inp);
              }
              return out;
          },
          "inp"_a);

    m.def("sum", [] (std::vector<std::tuple<double, FunctionTree<D> *>> &inp) {
              auto out = std::unique_ptr<FunctionTree<D>>(nullptr);
              if (inp.size() > 0) {
                  auto &mra = std::get<1>(inp[0])->getMRA();
                  out = std::make_unique<FunctionTree<D>>(mra);
                  build_grid(*out, inp);
                  add(-1.0, *out, inp);
              }
              return out;
          },
          "inp"_a);

    m.def("dot",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&mrcpp::dot<D>),
          "bra"_a,
          "ket"_a);

    m.def("prod", [] (std::vector<FunctionTree<D> *> &inp) {
              auto out = std::unique_ptr<FunctionTree<D>>(nullptr);
              if (inp.size() > 0) {
                  auto &mra = inp[0]->getMRA();
                  out = std::make_unique<FunctionTree<D>>(mra);
                  build_grid(*out, inp); // Union grid
                  build_grid(*out, 1);   // One extra refinement
                  multiply(-1.0, *out, inp);
              }
              return out;
          },
          "inp"_a);

    m.def("prod", [] (std::vector<std::tuple<double, FunctionTree<D> *>> &inp) {
              auto out = std::unique_ptr<FunctionTree<D>>(nullptr);
              if (inp.size() > 0) {
                  auto &mra = std::get<1>(inp[0])->getMRA();
                  out = std::make_unique<FunctionTree<D>>(mra);
                  build_grid(*out, inp); // Union grid
                  build_grid(*out, 1);   // One extra refinement
                  multiply(-1.0, *out, inp);
              }
              return out;
          },
          "inp"_a);
}

template <int D> void advanced_arithmetics(pybind11::module &m) {
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
