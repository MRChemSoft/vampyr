#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/treebuilders/apply.h>
//#include <MRCPP/treebuilders/complex_apply.h>

namespace vampyr {
template <int D> void applys(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "divergence",
        [](DerivativeOperator<D> &oper, std::vector<FunctionTree<D> *> &inp) {
            std::unique_ptr<FunctionTree<D>> out{nullptr};
            if (inp.size() == D) {
                out = std::make_unique<FunctionTree<D>>(inp[0]->getMRA());
                divergence<D>(*out, oper, inp);
            }
            return out;
        },
        "oper"_a,
        "inp"_a);

    m.def(
        "gradient",
        [](DerivativeOperator<D> &oper, FunctionTree<D> &inp) {
            auto tmp = mrcpp::gradient(oper, inp);
            std::vector<std::unique_ptr<FunctionTree<D>>> out;
            for (int i = 0; i < tmp.size(); i++) {
                auto *tmp_p = std::get<1>(tmp[i]);
                out.push_back(std::unique_ptr<FunctionTree<D>>(tmp_p));
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
        py::overload_cast<double, FunctionTree<D> &, ConvolutionOperator<D> &, FunctionTree<D> &, int, bool>(&apply<D>),
        "prec"_a,
        "out"_a,
        "oper"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def("apply",
          py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTree<D> &, int>(&apply<D>),
          "out"_a,
          "oper"_a,
          "inp"_a,
          "dir"_a = -1);
}

} // namespace vampyr
