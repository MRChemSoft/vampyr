#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/treebuilders/apply.h>

namespace vampyr {
template <int D> void applys(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "divergence",
        [](DerivativeOperator<D> &oper, std::vector<FunctionTree<D, double> *> &inp) {
            std::unique_ptr<FunctionTree<D, double>> out{nullptr};
            if (inp.size() == (size_t)D) {
                out = std::make_unique<FunctionTree<D, double>>(inp[0]->getMRA());
                divergence<D, double>(*out, oper, inp);
            }
            return out;
        },
        "oper"_a,
        "inp"_a);

    m.def(
        "gradient",
        [](DerivativeOperator<D> &oper, FunctionTree<D, double> &inp) {
            auto tmp = mrcpp::gradient<D, double>(oper, inp);
            std::vector<std::unique_ptr<FunctionTree<D, double>>> out;
            for (size_t i = 0; i < tmp.size(); i++) {
                auto *tmp_p = std::get<1>(tmp[i]);
                out.push_back(std::unique_ptr<FunctionTree<D, double>>(tmp_p));
            }
            mrcpp::clear(tmp, false);
            return out;
        },
        "oper"_a,
        "inp"_a);
}

// Direct bindings to MRCPP functionality
template <int D> void advanced_applys(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "apply",
        [](double prec, FunctionTree<D, double> &out, ConvolutionOperator<D> &oper, FunctionTree<D, double> &inp, int max_iter, bool abs_prec) {
            mrcpp::apply<D, double>(prec, out, oper, inp, max_iter, abs_prec);
        },
        "prec"_a,
        "out"_a,
        "oper"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def("apply",
          [](FunctionTree<D, double> &out, DerivativeOperator<D> &oper, FunctionTree<D, double> &inp, int dir) {
              mrcpp::apply<D, double>(out, oper, inp, dir);
          },
          "out"_a,
          "oper"_a,
          "inp"_a,
          "dir"_a = -1);
}

} // namespace vampyr
