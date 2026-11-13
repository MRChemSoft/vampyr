#pragma once

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

#include <string>

#include <MRCPP/operators/CartesianConvolution.h>
#include <MRCPP/operators/HelmholtzOperator.h>
#include <MRCPP/operators/IdentityConvolution.h>
#include <MRCPP/operators/PoissonOperator.h>
#include <MRCPP/operators/TimeEvolutionOperator.h>
#include <MRCPP/operators/HeatOperator.h>
#include <MRCPP/treebuilders/apply.h>

#include "trees/complex_utils.h"

namespace vampyr {

// MRCPP aborts on each of these, except a time <= 0, which gives a zero operator.
// It reads max_Jpower / 2 + 1 of its 60 cross-correlation matrices, hence 119.
template <int D>
void check_time_evolution_args(const mrcpp::MultiResolutionAnalysis<D> &mra,
                               double time,
                               int finest_scale,
                               int max_Jpower) {
    namespace py = pybind11;

    const auto &box = mra.getWorldBox();
    if (mra.getScalingBasis().getScalingType() != mrcpp::Legendre) {
        throw py::value_error("TimeEvolutionOperator requires a Legendre basis.");
    }
    if (box.isPeriodic()) throw py::value_error("TimeEvolutionOperator does not support a periodic world.");
    for (int d = 1; d < D; d++) {
        if (box.getScalingFactor(d) != box.getScalingFactor(0)) {
            throw py::value_error("TimeEvolutionOperator requires the same scaling factor in every direction.");
        }
    }
    if (!(time > 0.0)) throw py::value_error("time must be positive.");
    if (max_Jpower < 1 || max_Jpower > 119) {
        throw py::value_error("max_Jpower must be between 1 and 119, got " + std::to_string(max_Jpower) + ".");
    }
    if (finest_scale != mrcpp::TimeEvolutionOperator<D>::Adaptive && finest_scale < mra.getRootScale()) {
        throw py::value_error("finest_scale " + std::to_string(finest_scale) + " is coarser than the root scale " +
                              std::to_string(mra.getRootScale()) + ".");
    }
}

void cartesian_convolution(pybind11::module &);
void helmholtz_operator(pybind11::module &);
void poisson_operator(pybind11::module &);
template <int D> void time_evolution_operator(pybind11::module &m);
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
            "True if the operator has real coefficients")
        .def(
            "iscomplex",
            [](ConvolutionOperator<D> &C) { return static_cast<bool>(C.iscomplex()); },
            "True if the operator has complex coefficients")
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
    time_evolution_operator<D>(m);
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
        .def(
            "__call__",
            [](CartesianConvolution &O, FunctionTree<3, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<3, ComplexDouble>>(inp->getMRA());
                apply<3, ComplexDouble>(O.getBuildPrec(), *out, O, *inp);
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


template <int D> void time_evolution_operator(pybind11::module &m)
{
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    auto teo = py::class_<TimeEvolutionOperator<D>, ConvolutionOperator<D>>(m,
                                                                            "TimeEvolutionOperator",
                                                                            // clang-format off
    R"mydelimiter(
        The free-particle Schrodinger propagator exp(i t Laplacian), t > 0.
        For H = -Laplacian / 2 a time step dt is time = dt / 2.

        Applies to real and complex trees and returns a ComplexFunctionTree.
        The basis must be Legendre and the world non-periodic with the same
        scaling factor in every direction, for example BoundingBox([-30, 30]).

        With finest_scale the operator is built uniformly down to that scale;
        the default, TimeEvolutionOperator.Adaptive, builds it adaptively.
    )mydelimiter");
    // clang-format on

    if constexpr (D == 1) {
        // MRCPP deletes the signatures with the old 'imaginary' flag; here they
        // raise. Registered first: pybind11 takes a bool for an int, so the flag
        // would otherwise bind as finest_scale or max_Jpower.
        teo.def(py::init([](const MultiResolutionAnalysis<1> &, double, double, bool, int)
                             -> std::unique_ptr<TimeEvolutionOperator<1>> {
                    throw py::type_error("TimeEvolutionOperator no longer takes 'imaginary': the operator is complex. "
                                         "Use TimeEvolutionOperator(mra, prec, time).");
                }),
                "mra"_a,
                "prec"_a,
                "time"_a,
                "imaginary"_a,
                "max_Jpower"_a = 30,
                "Old signature with the removed 'imaginary' flag. Raises TypeError.");
        teo.def(py::init([](const MultiResolutionAnalysis<1> &, double, double, int, bool, int)
                             -> std::unique_ptr<TimeEvolutionOperator<1>> {
                    throw py::type_error("TimeEvolutionOperator no longer takes 'imaginary': the operator is complex. "
                                         "Use TimeEvolutionOperator(mra, prec, time, finest_scale).");
                }),
                "mra"_a,
                "prec"_a,
                "time"_a,
                "finest_scale"_a,
                "imaginary"_a,
                "max_Jpower"_a = 30,
                "Old signature with the removed 'imaginary' flag. Raises TypeError.");
    }

    teo.def(py::init(
                [](const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, int max_Jpower) {
                    check_time_evolution_args<D>(mra, time, finest_scale, max_Jpower);
                    return std::make_unique<TimeEvolutionOperator<D>>(mra, prec, time, finest_scale, max_Jpower);
                }),
            "mra"_a,
            "prec"_a,
            "time"_a,
            py::arg_v("finest_scale", TimeEvolutionOperator<D>::Adaptive, "TimeEvolutionOperator.Adaptive"),
            "max_Jpower"_a = 30)
        .def(
            "__call__",
            [](TimeEvolutionOperator<D> &T, FunctionTree<D, ComplexDouble> *inp) {
                auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
                apply<D, ComplexDouble>(T.getBuildPrec(), *out, T, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](TimeEvolutionOperator<D> &T, FunctionTree<D, double> *inp) {
                auto cinp = promote<D>(*inp);
                auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(cinp->getMRA());
                apply<D, ComplexDouble>(T.getBuildPrec(), *out, T, *cinp);
                return out;
            },
            "inp"_a,
            "Promotes the input; returns a ComplexFunctionTree");

    teo.attr("Adaptive") = TimeEvolutionOperator<D>::Adaptive;
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
