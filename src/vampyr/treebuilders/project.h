#pragma once

#include <pybind11/functional.h>
#include "PyProjectors.h"

namespace vampyr {
template <int D> void project(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("ZeroTree",
        [](const MultiResolutionAnalysis<D> &mra,
           const std::string &name) {
            auto out = std::make_unique<FunctionTree<D>>(mra, name);
            out->setZero();
            return out;
        },
        "mra"_a,
        "name"_a = "nn");

    py::class_<PyScalingProjector<D>>(m, "ScalingProjector")
        .def(py::init<const MultiResolutionAnalysis<D> &, double>(), "mra"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "scale"_a)
        .def(
            "__call__", [](PyScalingProjector<D> &P, RepresentableFunction<D> &func) { return P(func); }, "func"_a)
        .def(
            "__call__",
            [](PyScalingProjector<D> &P, std::function<double(const Coord<D> &r)> func) {

                try {
                    // When the analytic function func is badly defined, it kills the kernel
                    // of Notebooks. This evaluates func in a point, and if it is not successful
                    // it throws an error instead of killing the kernel.
                    auto arr = std::array<double, D>();
                    arr.fill(111111.111); // A number which hopefully does not divide by zero
                    func(arr);
                } catch (py::cast_error &e) {
                    py::print("Error: Invalid definition of analytic function");
                    throw;
                }
                auto old_threads = mrcpp_get_num_threads();
                set_max_threads(1);
                auto out = P(func);
                set_max_threads(old_threads);
                return out;
            },
            "func"_a);

    py::class_<PyWaveletProjector<D>>(m, "WaveletProjector")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "scale"_a)
        .def(
            "__call__", [](PyWaveletProjector<D> &P, RepresentableFunction<D> &func) { return P(func); }, "func"_a)
        .def(
            "__call__",
            [](PyWaveletProjector<D> &P, std::function<double(const Coord<D> &r)> func) {

                try {
                    auto arr = std::array<double, D>();
                    arr.fill(111111.111); // A number which hopefully does not divide by zero
                    func(arr);
                } catch (py::cast_error &e) {
                    py::print("Error: Invalid definition of analytic function");
                    throw;
                }

                auto old_threads = mrcpp_get_num_threads();
                set_max_threads(1);
                auto out = P(func);
                set_max_threads(old_threads);
                return out;
            },
            "func"_a);
}

template <int D> void advanced_project(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("project",
          py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int, bool>(&mrcpp::project<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def(
        "project",
        [](double prec,
           FunctionTree<D> &out,
           std::function<double(const Coord<D> &r)> inp,
           int max_iter,
           bool abs_prec) {
            auto old_threads = mrcpp_get_num_threads();
            mrcpp::set_max_threads(1);
            mrcpp::project<D>(prec, out, inp, max_iter, abs_prec);
            mrcpp::set_max_threads(old_threads);
        },
        "prec"_a = -1.0,
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);
}
} // namespace vampyr
