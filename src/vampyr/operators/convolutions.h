#pragma once

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

#include <MRCPP/operators/CartesianConvolution.h>
#include <MRCPP/operators/HelmholtzOperator.h>
#include <MRCPP/operators/IdentityConvolution.h>
#include <MRCPP/operators/PoissonOperator.h>
#include <MRCPP/operators/TimeEvolutionOperator.h>
#include <MRCPP/operators/HeatOperator.h>
#include <MRCPP/treebuilders/apply.h>

#include <string>

#include "../trees/complex_utils.h"

namespace vampyr {

// The power integrals are built for scales 0..N+1, and the operator root is
// the MRA root scale, so a negative root scale dereferences a null map entry
// in the cross correlation calculator. The time is scaled as time * 4^n with
// no scaling factor, so a stretched world evolves for the wrong duration.
template <int D> void check_time_evolution_mra(const mrcpp::MultiResolutionAnalysis<D> &mra) {
    namespace py = pybind11;

    if (mra.getRootScale() < 0) {
        throw py::value_error("TimeEvolutionOperator requires a non-negative root scale, got " +
                              std::to_string(mra.getRootScale()) + ". Use BoundingBox(0).");
    }
    for (int d = 0; d < D; d++) {
        if (mra.getWorldBox().getScalingFactor(d) != 1.0) {
            throw py::value_error("TimeEvolutionOperator does not account for the world scaling "
                                  "factor. Use a world without one.");
        }
    }
}


void cartesian_convolution(pybind11::module &);
void helmholtz_operator(pybind11::module &);
void poisson_operator(pybind11::module &);
void time_evolution_operator(pybind11::module &m);
void heat_operator(pybind11::module &m);

template <int D> void convolutions(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<ConvolutionOperator<D>>(m, "ConvolutionOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &, GaussExp<1> &, double>(), "mra"_a, "kernel"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<D> &, GaussExp<1> &, double, int, int>())
        .def(
            "isreal",
            [](ConvolutionOperator<D> &C) { return static_cast<bool>(C.isreal()); },
            "Real coefficients")
        .def(
            "iscomplex",
            [](ConvolutionOperator<D> &C) { return static_cast<bool>(C.iscomplex()); },
            "Complex coefficients; will not apply to a real FunctionTree")
        .def(
            "__call__",
            [](ConvolutionOperator<D> &C, FunctionTree<D, double> *inp) {
                auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
                apply<D, double>(C.getBuildPrec(), *out, C, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](ConvolutionOperator<D> &C, FunctionTree<D, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
                apply<D, ComplexDouble>(C.getBuildPrec(), *out, C, *inp);
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
            [](IdentityConvolution<D> &I, FunctionTree<D, double> *inp) {
                auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
                apply<D, double>(I.getBuildPrec(), *out, I, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](IdentityConvolution<D> &I, FunctionTree<D, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
                apply<D, ComplexDouble>(I.getBuildPrec(), *out, I, *inp);
                return out;
            },
            "inp"_a);

    if constexpr (D == 3) cartesian_convolution(m);
    if constexpr (D == 3) helmholtz_operator(m);
    if constexpr (D == 3) poisson_operator(m);
    if constexpr (D == 1) time_evolution_operator(m);
    if constexpr (D == 1) heat_operator(m);
}

void cartesian_convolution(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<CartesianConvolution, ConvolutionOperator<3>>(m, "CartesianConvolution")
        .def(py::init<const MultiResolutionAnalysis<3> &, GaussExp<1> &, double>(), "mra"_a, "kernel"_a, "prec"_a)
        .def(
            "__call__",
            [](CartesianConvolution &O, FunctionTree<3, double> *inp) {
                auto out = std::make_unique<FunctionTree<3, double>>(inp->getMRA());
                apply<3, double>(O.getBuildPrec(), *out, O, *inp);
                return out;
            },
            "inp"_a)
        .def("setCartesianComponents", &CartesianConvolution::setCartesianComponents);
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
            [](PoissonOperator &P, FunctionTree<3, double> *inp) {
                auto out = std::make_unique<FunctionTree<3, double>>(inp->getMRA());
                apply<3, double>(P.getBuildPrec(), *out, P, *inp);
                out->rescale(1.0 / (4.0 * mrcpp::pi));
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](PoissonOperator &P, FunctionTree<3, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<3, ComplexDouble>>(inp->getMRA());
                apply<3, ComplexDouble>(P.getBuildPrec(), *out, P, *inp);
                out->rescale(ComplexDouble(1.0 / (4.0 * mrcpp::pi), 0.0));
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
            [](HelmholtzOperator &H, FunctionTree<3, double> *inp) {
                auto out = std::make_unique<FunctionTree<3, double>>(inp->getMRA());
                apply<3, double>(H.getBuildPrec(), *out, H, *inp);
                out->rescale(1.0 / (4.0 * mrcpp::pi));
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](HelmholtzOperator &H, FunctionTree<3, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<3, ComplexDouble>>(inp->getMRA());
                apply<3, ComplexDouble>(H.getBuildPrec(), *out, H, *inp);
                out->rescale(ComplexDouble(1.0 / (4.0 * mrcpp::pi), 0.0));
                return out;
            },
            "inp"_a);
}


void time_evolution_operator(pybind11::module &m)
{
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    auto teo = py::class_<TimeEvolutionOperator<1>, ConvolutionOperator<1>>(m,
                                                                "TimeEvolutionOperator",
                                                                // clang-format off
    R"mydelimiter(
        The free particle Schrodinger semigroup exp(i t d^2/dx^2).

        One operator tree carries the whole cos + 1j*sin kernel, so there is no
        imaginary flag and no 2x2 real block composition. Returns a
        ComplexFunctionTree; a real input is promoted on the way in.

        finest_scale refines uniformly down to that scale; the default,
        TimeEvolutionOperator.Adaptive, refines adaptively.
    )mydelimiter")
        // clang-format on
        // The two rejected legacy forms are registered ahead of the current
        // fixed-scale overload on purpose. Python bools satisfy pybind11's
        // integer caster even in the no-convert pass, so without these an old
        // call like U(mra, prec, t, scale, False) would bind silently with
        // max_Jpower = 0 instead of failing.
        .def(py::init([](const MultiResolutionAnalysis<1> &, double, double, bool, int)
                          -> std::unique_ptr<TimeEvolutionOperator<1>> {
                 throw py::type_error("The 'imaginary' flag is gone, the kernel carries both parts. "
                                      "Use TimeEvolutionOperator(mra, prec, time).");
             }),
             "mra"_a,
             "prec"_a,
             "time"_a,
             "imaginary"_a,
             "max_Jpower"_a = 30)
        .def(py::init([](const MultiResolutionAnalysis<1> &, double, double, int, bool, int)
                          -> std::unique_ptr<TimeEvolutionOperator<1>> {
                 throw py::type_error("The 'imaginary' flag is gone, the kernel carries both parts. "
                                      "Use TimeEvolutionOperator(mra, prec, time, finest_scale).");
             }),
             "mra"_a,
             "prec"_a,
             "time"_a,
             "finest_scale"_a,
             "imaginary"_a,
             "max_Jpower"_a = 30)
        .def(py::init([](const MultiResolutionAnalysis<1> &mra,
                         double prec,
                         double time,
                         int finest_scale,
                         int max_Jpower) {
                 check_time_evolution_mra<1>(mra);
                 return std::make_unique<TimeEvolutionOperator<1>>(mra, prec, time, finest_scale, max_Jpower);
             }),
             "mra"_a,
             "prec"_a,
             "time"_a,
             "finest_scale"_a = TimeEvolutionOperator<1>::Adaptive,
             "max_Jpower"_a = 30)
        .def(
            "__call__",
            [](TimeEvolutionOperator<1> &T, FunctionTree<1, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<1, ComplexDouble>>(inp->getMRA());
                apply<1, ComplexDouble>(T.getBuildPrec(), *out, T, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](TimeEvolutionOperator<1> &T, FunctionTree<1, double> *inp) {
                auto cinp = promote<1>(*inp);
                auto out = std::make_unique<FunctionTree<1, ComplexDouble>>(cinp->getMRA());
                apply<1, ComplexDouble>(T.getBuildPrec(), *out, T, *cinp);
                return out;
            },
            "inp"_a,
            "Promotes the input; returns a ComplexFunctionTree");

    teo.attr("Adaptive") = TimeEvolutionOperator<1>::Adaptive;
}


void heat_operator(pybind11::module &m)
{
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    py::class_<HeatOperator<1>, ConvolutionOperator<1>>(m, "HeatOperator")
        .def(py::init<const MultiResolutionAnalysis<1> &, double, double>(),
             "mra"_a,
             "time"_a,
             "prec"_a)
        .def(
            "__call__",
            [](HeatOperator<1> &T, FunctionTree<1, double> *inp) {
                auto out = std::make_unique<FunctionTree<1, double>>(inp->getMRA());
                apply<1, double>(T.getBuildPrec(), *out, T, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](HeatOperator<1> &T, FunctionTree<1, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<1, ComplexDouble>>(inp->getMRA());
                apply<1, ComplexDouble>(T.getBuildPrec(), *out, T, *inp);
                return out;
            },
            "inp"_a);
}

} // namespace vampyr
