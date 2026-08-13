#pragma once

/*
 * Complex-valued bindings for VAMPyR.
 *
 * Exposes the T = ComplexDouble instantiations that MRCPP provides since the
 * v2 templatization of FunctionTree<D, T> and the treebuilders. Everything in
 * this header is additive: real-valued bindings are untouched, mixed
 * real/complex arithmetic is resolved through Python's reflected operator
 * protocol (__radd__ & co.), relying on py::is_operator() returning
 * NotImplemented on overload misses.
 */

#include <complex>
#include <filesystem>
#include <functional>
#include <memory>

#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <MRCPP/operators/TimeEvolutionOperator.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWTree.h>
#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/apply.h>
#include <MRCPP/treebuilders/grid.h>
#include <MRCPP/treebuilders/multiply.h>
#include <MRCPP/treebuilders/project.h>
#include <MRCPP/treebuilders/treeMap.h>

#include "../treebuilders/PyFunctionMap.h"
#include "../treebuilders/PyProjectors.h"

namespace vampyr {

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

/** Deep copy of a real tree into a complex tree (exact, same grid). */
template <int D>
auto promote(mrcpp::FunctionTree<D, double> &inp) -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    return std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>>(inp.CopyTreeToComplex());
}

template <int D>
auto cplx_linear_comb(ComplexDouble a,
                      mrcpp::FunctionTree<D, ComplexDouble> &inp_a,
                      ComplexDouble b,
                      mrcpp::FunctionTree<D, ComplexDouble> &inp_b)
    -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp_a.getMRA());
    FunctionTreeVector<D, ComplexDouble> vec;
    vec.push_back({a, &inp_a});
    vec.push_back({b, &inp_b});
    build_grid(*out, vec);
    add(-1.0, *out, vec);
    return out;
}

template <int D>
auto cplx_scale(ComplexDouble c, mrcpp::FunctionTree<D, ComplexDouble> &inp)
    -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp.getMRA());
    copy_grid(*out, inp);
    copy_func(*out, inp);
    out->rescale(c);
    return out;
}

template <int D>
auto cplx_multiply(mrcpp::FunctionTree<D, ComplexDouble> &inp_a, mrcpp::FunctionTree<D, ComplexDouble> &inp_b)
    -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp_a.getMRA());
    FunctionTreeVector<D, ComplexDouble> vec;
    vec.push_back({ComplexDouble(1.0, 0.0), &inp_a});
    vec.push_back({ComplexDouble(1.0, 0.0), &inp_b});
    build_grid(*out, vec);
    build_grid(*out, 1);
    multiply(-1.0, *out, vec);
    return out;
}

/** Materialized complex conjugate (exact, same grid). */
template <int D>
auto cplx_conj(mrcpp::FunctionTree<D, ComplexDouble> &inp) -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp.getMRA());
    copy_grid(*out, inp);
    treeMap<D, ComplexDouble>(-1.0, *out, inp, [](ComplexDouble z) { return std::conj(z); });
    return out;
}

/* ------------------------------------------------------------------ */
/* Trees                                                              */
/* ------------------------------------------------------------------ */

template <int D> void trees_complex(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<MWTree<D, ComplexDouble>>(m,
                                        "ComplexMWTree",
                                        // clang-format off
    R"mydelimiter(
        Base class for complex-valued multiwavelet trees. Same interface as
        :class:`MWTree`, with complex scaling and wavelet coefficients.
    )mydelimiter")
        // clang-format on
        .def("MRA", &MWTree<D, ComplexDouble>::getMRA, py::return_value_policy::reference_internal)
        .def("nNodes", &MWTree<D, ComplexDouble>::getNNodes)
        .def("nEndNodes", &MWTree<D, ComplexDouble>::getNEndNodes)
        .def("nRootNodes", &MWTree<D, ComplexDouble>::getNRootNodes)
        .def("rootScale", &MWTree<D, ComplexDouble>::getRootScale)
        .def("depth", &MWTree<D, ComplexDouble>::getDepth)
        .def("setZero",
             [](MWTree<D, ComplexDouble> *out) {
                 out->setZero();
                 return out;
             })
        .def("clear", &MWTree<D, ComplexDouble>::clear)
        .def("setName", &MWTree<D, ComplexDouble>::setName)
        .def("name", &MWTree<D, ComplexDouble>::getName)
        .def("squaredNorm", &MWTree<D, ComplexDouble>::getSquareNorm)
        .def("norm", [](MWTree<D, ComplexDouble> &tree) {
            auto sqNorm = tree.getSquareNorm();
            return (sqNorm >= 0.0) ? std::sqrt(sqNorm) : -1.0;
        });

    py::class_<FunctionTree<D, ComplexDouble>, MWTree<D, ComplexDouble>>(m,
                                                                        "ComplexFunctionTree",
                                                                        // clang-format off
    R"mydelimiter(
        A complex-valued function in the multiwavelet basis.

        Arithmetic works as for the real :class:`FunctionTree`, and mixes
        freely with real trees and Python scalars: a real operand is promoted
        on the spot. Use :meth:`real`, :meth:`imag` and :meth:`conj` to move
        between the two, and ``from_real_imag`` to build one from a pair of
        real trees.
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, const std::string &>(), "mra"_a, "name"_a = "nn")
        .def(py::init([](FunctionTree<D, double> &re, FunctionTree<D, double> &im) {
                 return std::make_unique<FunctionTree<D, ComplexDouble>>(re, im);
             }),
             "real"_a,
             "imag"_a,
             "Construct the complex function real + 1j*imag from two real trees")
        .def("integrate", &FunctionTree<D, ComplexDouble>::integrate)
        .def("normalize",
             [](FunctionTree<D, ComplexDouble> *out) {
                 out->normalize();
                 return out;
             })
        .def("rescale", [](FunctionTree<D, ComplexDouble> *out, ComplexDouble c) { out->rescale(c); })
        .def(
            "crop",
            [](FunctionTree<D, ComplexDouble> *out, double prec, bool abs_prec) {
                out->crop(prec, 1.0, abs_prec);
                return out;
            },
            "prec"_a,
            "abs_prec"_a = false)
        .def("deepCopy",
             [](FunctionTree<D, ComplexDouble> *inp) {
                 auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
                 copy_grid(*out, *inp);
                 copy_func(*out, *inp);
                 return out;
             })
        .def("real",
             [](FunctionTree<D, ComplexDouble> &inp) { return std::unique_ptr<FunctionTree<D, double>>(inp.Real()); },
             "Real part as a new real FunctionTree")
        .def("imag",
             [](FunctionTree<D, ComplexDouble> &inp) { return std::unique_ptr<FunctionTree<D, double>>(inp.Imag()); },
             "Imaginary part as a new real FunctionTree")
        .def("conj", [](FunctionTree<D, ComplexDouble> &inp) { return cplx_conj<D>(inp); },
             "Complex conjugate as a new ComplexFunctionTree (exact, same grid)")
        .def(
            "saveTree",
            [](FunctionTree<D, ComplexDouble> &obj, const std::string &filename) {
                namespace fs = std::filesystem;
                obj.saveTree(filename);
                return fs::absolute(fs::path(filename + ".tree"));
            },
            "filename"_a)
        .def("loadTree", &FunctionTree<D, ComplexDouble>::loadTree, "filename"_a)
        .def("__call__",
             [](FunctionTree<D, ComplexDouble> &func, const Coord<D> &r) { return func.evalf_precise(r); })
        .def("__pos__",
             [](FunctionTree<D, ComplexDouble> *inp) {
                 auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
                 copy_grid(*out, *inp);
                 copy_func(*out, *inp);
                 return out;
             },
             py::is_operator())
        .def("__neg__",
             [](FunctionTree<D, ComplexDouble> *inp) { return cplx_scale<D>(ComplexDouble(-1.0, 0.0), *inp); },
             py::is_operator())
        // complex (+|-) complex
        .def("__add__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) {
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {1.0, 0.0}, *b);
             },
             py::is_operator())
        .def("__iadd__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) {
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {1.0, 0.0}, *b);
             },
             py::is_operator())
        .def("__sub__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) {
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {-1.0, 0.0}, *b);
             },
             py::is_operator())
        .def("__isub__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) {
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {-1.0, 0.0}, *b);
             },
             py::is_operator())
        // complex (+|-) real (and reflected forms, hit when the real tree's
        // own operators return NotImplemented)
        .def("__add__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {1.0, 0.0}, *bc);
             },
             py::is_operator())
        .def("__radd__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_linear_comb<D>({1.0, 0.0}, *bc, {1.0, 0.0}, *a);
             },
             py::is_operator())
        .def("__sub__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_linear_comb<D>({1.0, 0.0}, *a, {-1.0, 0.0}, *bc);
             },
             py::is_operator())
        .def("__rsub__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_linear_comb<D>({1.0, 0.0}, *bc, {-1.0, 0.0}, *a);
             },
             py::is_operator())
        // products
        .def("__mul__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) { return cplx_multiply<D>(*a, *b); },
             py::is_operator())
        .def("__imul__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, ComplexDouble> *b) { return cplx_multiply<D>(*a, *b); },
             py::is_operator())
        .def("__mul__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_multiply<D>(*a, *bc);
             },
             py::is_operator())
        .def("__rmul__",
             [](FunctionTree<D, ComplexDouble> *a, FunctionTree<D, double> *b) {
                 auto bc = promote<D>(*b);
                 return cplx_multiply<D>(*bc, *a);
             },
             py::is_operator())
        .def("__mul__",
             [](FunctionTree<D, ComplexDouble> *a, ComplexDouble c) { return cplx_scale<D>(c, *a); },
             py::is_operator())
        .def("__imul__",
             [](FunctionTree<D, ComplexDouble> *a, ComplexDouble c) { return cplx_scale<D>(c, *a); },
             py::is_operator())
        .def("__rmul__",
             [](FunctionTree<D, ComplexDouble> *a, ComplexDouble c) { return cplx_scale<D>(c, *a); },
             py::is_operator())
        .def("__truediv__",
             [](FunctionTree<D, ComplexDouble> *a, ComplexDouble c) { return cplx_scale<D>(1.0 / c, *a); },
             py::is_operator())
        .def("__itruediv__",
             [](FunctionTree<D, ComplexDouble> *a, ComplexDouble c) { return cplx_scale<D>(1.0 / c, *a); },
             py::is_operator());

    m.def(
        "from_real_imag",
        [](FunctionTree<D, double> &re, FunctionTree<D, double> &im) {
            return std::make_unique<FunctionTree<D, ComplexDouble>>(re, im);
        },
        "real"_a,
        "imag"_a,
        "Build the complex function real + 1j*imag from two real FunctionTrees");
}

/* ------------------------------------------------------------------ */
/* Projection                                                         */
/* ------------------------------------------------------------------ */

template <int D> void project_complex(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<PyScalingProjector<D, ComplexDouble>>(m,
                                                     "ComplexScalingProjector",
                                                     // clang-format off
    R"mydelimiter(
        Projects a complex-valued Python callable onto the scaling basis,
        returning a :class:`ComplexFunctionTree`. Normally reached as
        ``ScalingProjector(mra, prec, dtype=complex)``.
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, double>(), "mra"_a, "prec"_a)
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "scale"_a)
        .def(
            "__call__",
            [](PyScalingProjector<D, ComplexDouble> &P, std::function<ComplexDouble(const Coord<D> &r)> func) {
                try {
                    auto arr = std::array<double, D>();
                    arr.fill(111111.111);
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

/* ------------------------------------------------------------------ */
/* Arithmetic free functions + map                                    */
/* ------------------------------------------------------------------ */

template <int D> void arithmetics_complex(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // NOTE: for complex trees mrcpp::dot conjugates the bra by default, i.e.
    // dot(bra, ket) is the L2 inner product <bra|ket>. The plain bilinear
    // integral of a product is dot(bra.conj(), ket).
    m.def(
        "dot",
        [](FunctionTree<D, ComplexDouble> &bra, FunctionTree<D, ComplexDouble> &ket) {
            return mrcpp::dot<D, ComplexDouble, ComplexDouble>(bra, ket);
        },
        "bra"_a,
        "ket"_a,
        "L2 inner product <bra|ket>; the bra is conjugated. For the bilinear "
        "integral of the product use dot(bra.conj(), ket)");
    m.def(
        "dot",
        [](FunctionTree<D, ComplexDouble> &bra, FunctionTree<D, double> &ket) {
            return mrcpp::dot<D, ComplexDouble, double>(bra, ket);
        },
        "bra"_a,
        "ket"_a);
    m.def(
        "dot",
        [](FunctionTree<D, double> &bra, FunctionTree<D, ComplexDouble> &ket) {
            return mrcpp::dot<D, double, ComplexDouble>(bra, ket);
        },
        "bra"_a,
        "ket"_a);

    py::class_<PyFunctionMap<D, ComplexDouble>>(m,
                                                "ComplexFunctionMap",
                                                // clang-format off
    R"mydelimiter(
        Applies a pointwise complex map through the multiwavelet
        representation. Normally reached as
        ``FunctionMap(fmap, prec, dtype=complex)``.
    )mydelimiter")
        // clang-format on
        .def(py::init<std::function<ComplexDouble(ComplexDouble)>, double>(), "fmap"_a, "prec"_a)
        .def(
            "__call__",
            [](PyFunctionMap<D, ComplexDouble> &F, FunctionTree<D, ComplexDouble> &inp) {
                auto old_threads = mrcpp_get_num_threads();
                set_max_threads(1);
                auto out = F(inp);
                set_max_threads(old_threads);
                return out;
            },
            "inp"_a);
}

/* ------------------------------------------------------------------ */
/* advanced.* overloads                                               */
/* ------------------------------------------------------------------ */

template <int D> void advanced_complex(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "apply",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           ConvolutionOperator<D> &oper,
           FunctionTree<D, ComplexDouble> &inp,
           int max_iter,
           bool abs_prec) { mrcpp::apply<D, ComplexDouble>(prec, out, oper, inp, max_iter, abs_prec); },
        "prec"_a,
        "out"_a,
        "oper"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def(
        "apply",
        [](FunctionTree<D, ComplexDouble> &out, DerivativeOperator<D> &oper, FunctionTree<D, ComplexDouble> &inp, int dir) {
            mrcpp::apply<D, ComplexDouble>(out, oper, inp, dir);
        },
        "out"_a,
        "oper"_a,
        "inp"_a,
        "dir"_a = -1);

    m.def(
        "add",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           ComplexDouble a,
           FunctionTree<D, ComplexDouble> &inp_a,
           ComplexDouble b,
           FunctionTree<D, ComplexDouble> &inp_b,
           int max_iter,
           bool abs_prec) { mrcpp::add<D, ComplexDouble>(prec, out, a, inp_a, b, inp_b, max_iter, abs_prec); },
        "prec"_a = -1.0,
        "out"_a,
        "a"_a = 1.0,
        "inp_a"_a,
        "b"_a = 1.0,
        "inp_b"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def(
        "multiply",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           ComplexDouble c,
           FunctionTree<D, ComplexDouble> &inp_a,
           FunctionTree<D, ComplexDouble> &inp_b,
           int max_iter,
           bool abs_prec,
           bool use_max_norms,
           bool conjugate) {
            mrcpp::multiply<D, ComplexDouble>(prec, out, c, inp_a, inp_b, max_iter, abs_prec, use_max_norms, conjugate);
        },
        "prec"_a = -1.0,
        "out"_a,
        "c"_a = 1.0,
        "inp_a"_a,
        "inp_b"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false,
        "use_max_norms"_a = false,
        "conjugate"_a = false);

    m.def(
        "project",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           std::function<ComplexDouble(const Coord<D> &r)> inp,
           int max_iter,
           bool abs_prec) {
            auto old_threads = mrcpp_get_num_threads();
            mrcpp::set_max_threads(1);
            mrcpp::project<D, ComplexDouble>(prec, out, inp, max_iter, abs_prec);
            mrcpp::set_max_threads(old_threads);
        },
        "prec"_a = -1.0,
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def(
        "map",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           FunctionTree<D, ComplexDouble> &inp,
           std::function<ComplexDouble(ComplexDouble)> fmap,
           int max_iter,
           bool abs_prec) {
            auto old_threads = mrcpp_get_num_threads();
            mrcpp::set_max_threads(1);
            mrcpp::treeMap<D, ComplexDouble>(prec, out, inp, fmap, max_iter, abs_prec);
            mrcpp::set_max_threads(old_threads);
        },
        "prec"_a = -1.0,
        "out"_a,
        "inp"_a,
        "fmap"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false);

    m.def(
        "build_grid",
        [](FunctionTree<D, ComplexDouble> &out, FunctionTree<D, ComplexDouble> &inp, int max_iter) {
            mrcpp::build_grid<D, ComplexDouble>(out, inp, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "copy_grid",
        [](FunctionTree<D, ComplexDouble> &out, FunctionTree<D, ComplexDouble> &inp) {
            mrcpp::copy_grid<D, ComplexDouble>(out, inp);
        },
        "out"_a,
        "inp"_a);

    m.def(
        "copy_func",
        [](FunctionTree<D, ComplexDouble> &out, FunctionTree<D, ComplexDouble> &inp) {
            mrcpp::copy_func<D, ComplexDouble>(out, inp);
        },
        "out"_a,
        "inp"_a);
}

/* ------------------------------------------------------------------ */
/* Complex Schrodinger time evolution (1D)                            */
/* ------------------------------------------------------------------ */

/** exp(i*t*d^2/dx^2) as a single complex-valued operator.
 *
 *  MRCPP builds one TimeEvolutionOperator<1, ComplexDouble> whose kernel
 *  carries the full cos + i*sin coefficient, so the semigroup is applied in
 *  one convolution rather than composed from two real operator trees. Only
 *  the fixed finest_scale constructor is instantiated for ComplexDouble in
 *  MRCPP, so the adaptive overload is deliberately not exposed here.
 *
 *  The `imaginary` flag of the real operator has no meaning for the complex
 *  kernel (the if-constexpr branch in TimeEvolution_CrossCorrelationCalculator
 *  emits both parts); it is pinned to false and kept off the Python API. */
inline void complex_time_evolution(pybind11::module &m) {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using ComplexTimeEvolution = mrcpp::TimeEvolutionOperator<1, ComplexDouble>;

    py::class_<ComplexTimeEvolution>(m,
                                     "ComplexTimeEvolutionOperator",
                                     // clang-format off
    R"mydelimiter(
        The free-particle Schrodinger semigroup exp(i t d^2/dx^2) as a single
        complex-valued convolution.

        The kernel carries the full cos + i*sin coefficient, so one
        application does what previously required two real operators composed
        as Re[U] psi + i Im[U] psi. Built at a fixed finest scale; there is no
        adaptive constructor.
    )mydelimiter")
        // clang-format on
        .def(py::init([](const mrcpp::MultiResolutionAnalysis<1> &mra,
                         double prec,
                         double time,
                         int finest_scale,
                         int max_Jpower) {
                 return std::make_unique<ComplexTimeEvolution>(mra, prec, time, finest_scale, false, max_Jpower);
             }),
             "mra"_a,
             "prec"_a,
             "time"_a,
             "finest_scale"_a,
             "max_Jpower"_a = 20)
        .def("buildPrec", &ComplexTimeEvolution::getBuildPrec)
        .def(
            "__call__",
            [](ComplexTimeEvolution &U, mrcpp::FunctionTree<1, ComplexDouble> *inp) {
                auto out = std::make_unique<mrcpp::FunctionTree<1, ComplexDouble>>(inp->getMRA());
                mrcpp::apply<1, ComplexDouble, ComplexDouble>(U.getBuildPrec(), *out, U, *inp);
                return out;
            },
            "inp"_a)
        .def(
            "__call__",
            [](ComplexTimeEvolution &U, mrcpp::FunctionTree<1, double> *inp) {
                auto cinp = promote<1>(*inp);
                auto out = std::make_unique<mrcpp::FunctionTree<1, ComplexDouble>>(cinp->getMRA());
                mrcpp::apply<1, ComplexDouble, ComplexDouble>(U.getBuildPrec(), *out, U, *cinp);
                return out;
            },
            "inp"_a);
}

/* ------------------------------------------------------------------ */
/* Entry point                                                        */
/* ------------------------------------------------------------------ */

template <int D> void bind_complex(pybind11::module &m) {
    namespace py = pybind11;

    trees_complex<D>(m);
    project_complex<D>(m);
    arithmetics_complex<D>(m);

    auto advanced = m.attr("advanced").cast<py::module>();
    advanced_complex<D>(advanced);

    if constexpr (D == 1) complex_time_evolution(m);
}

} // namespace vampyr
