#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/HelmholtzOperator.h>
#include <MRCPP/operators/IdentityConvolution.h>
#include <MRCPP/operators/PoissonOperator.h>
#include <MRCPP/treebuilders/apply.h>

namespace vampyr {

void helmholtz_operator(pybind11::module &);
void poisson_operator(pybind11::module &);

template <int D> void convolutions(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<ConvolutionOperator<D>>(m, "ConvolutionOperator")
    .def(py::init<const MultiResolutionAnalysis<D> &, GaussExp<1> &, double>())
    .def(py::init<const MultiResolutionAnalysis<D> &, GaussExp<1> &, double, int, int>())
    .def(
        "__call__",
        [](ConvolutionOperator<D> &C, FunctionTree<D> *inp) {
            auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
            apply<D>(C.getBuildPrec(), *out, C, *inp);
            return out;
            },
            "inp"_a);

    py::class_<IdentityConvolution<D>, ConvolutionOperator<D>>(m, "IdentityConvolution")
        .def(py::init<const MultiResolutionAnalysis<D> &, double>(), "mra"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<D> &, double, int, int>(),
             "mra"_a,
             "prec"_a,
             "root"_a = 0,
             "reach"_a = 1)
        .def(
            "__call__",
            [](IdentityConvolution<D> &I, FunctionTree<D> *inp) {
                auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
                apply<D>(I.getBuildPrec(), *out, I, *inp);
                return out;
            },
            "inp"_a);

    if constexpr (D == 3) helmholtz_operator(m);
    if constexpr (D == 3) poisson_operator(m);
}

void poisson_operator(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<PoissonOperator, ConvolutionOperator<3>>(m, "PoissonOperator")
        .def(py::init<const MultiResolutionAnalysis<3> &, double>(), "mra"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<3> &, double, int, int>(),
             "mra"_a,
             "prec"_a,
             "root"_a = 0,
             "reach"_a = 1)
        .def(
            "__call__",
            [](PoissonOperator &P, FunctionTree<3> *inp) {
                auto out = std::make_unique<FunctionTree<3>>(inp->getMRA());
                apply<3>(P.getBuildPrec(), *out, P, *inp);
                out->rescale(1.0 / (4.0 * mrcpp::pi));
                return out;
            },
            "inp"_a);
}

void helmholtz_operator(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<HelmholtzOperator, ConvolutionOperator<3>>(m, "HelmholtzOperator")
        .def(py::init<const MultiResolutionAnalysis<3> &, double, double>(), "mra"_a, "exp"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<3> &, double, double, int, int>(),
             "mra"_a,
             "exp"_a,
             "prec"_a,
             "root"_a = 0,
             "reach"_a = 1)
        .def(
            "__call__",
            [](HelmholtzOperator &H, FunctionTree<3> *inp) {
                auto out = std::make_unique<FunctionTree<3>>(inp->getMRA());
                apply<3>(H.getBuildPrec(), *out, H, *inp);
                out->rescale(1.0 / (4.0 * mrcpp::pi));
                return out;
            },
            "inp"_a);
}

} // namespace vampyr
