#pragma once

/*
 * Complex-valued bindings for VAMPyR.
 *
 * Exposes the T = ComplexDouble instantiations of FunctionTree<D, T> and the
 * treebuilders. Additive: the real bindings are untouched. Mixed real/complex
 * arithmetic goes through Python's reflected operators (__radd__ and friends),
 * which works because py::is_operator() returns NotImplemented on a miss.
 *
 * Operators are not bound here. MRCPP dispatches the operator scalar at
 * runtime, so the complex overloads sit next to their real counterparts in
 * operators/convolutions.h and operators/derivatives.h.
 */

#include <complex>
#include <sstream>
#include <filesystem>
#include <functional>
#include <memory>

#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

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
#include "../trees/complex_utils.h"

namespace vampyr {

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

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

/** Mirrors impl__pow__ in trees/trees.h. */
template <int D>
auto impl__pow__complex(mrcpp::FunctionTree<D, ComplexDouble> *inp, double c)
    -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    using namespace mrcpp;
    auto out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp->getMRA());
    copy_grid(*out, *inp);
    copy_func(*out, *inp);
    refine_grid(*out, 1);
    out->power(c);
    return out;
}

template <int D> void trees_complex(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<MWTree<D, ComplexDouble>>(m,
                                        "ComplexMWTree",
                                        // clang-format off
    R"mydelimiter(
        MWTree with complex scaling and wavelet coefficients.
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
        })
        .def("__str__",
             [](MWTree<D, ComplexDouble> &tree) {
                 std::ostringstream os;
                 os << tree;
                 return os.str();
             });

    py::class_<FunctionTree<D, ComplexDouble>, MWTree<D, ComplexDouble>>(m,
                                                                        "ComplexFunctionTree",
                                                                        // clang-format off
    R"mydelimiter(
        FunctionTree with complex coefficients.

        Arithmetic mixes with real trees and Python scalars, promoting the real
        operand. real(), imag() and conj() go back the other way.
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, const std::string &>(), "mra"_a, "name"_a = "nn")
        .def(py::init([](FunctionTree<D, double> &re, FunctionTree<D, double> &im) {
                 return std::make_unique<FunctionTree<D, ComplexDouble>>(re, im);
             }),
             "real"_a,
             "imag"_a,
             "Builds real + 1j*imag from two real trees")
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
             "Real part, as a new FunctionTree")
        .def("imag",
             [](FunctionTree<D, ComplexDouble> &inp) { return std::unique_ptr<FunctionTree<D, double>>(inp.Imag()); },
             "Imaginary part, as a new FunctionTree")
        .def("__pow__", &impl__pow__complex<D>, py::is_operator())
        .def("__ipow__", &impl__pow__complex<D>, py::is_operator())
        .def("conj", [](FunctionTree<D, ComplexDouble> &inp) { return cplx_conj<D>(inp); },
             "Complex conjugate, same grid")
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
        "Builds real + 1j*imag from two real trees");
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
        Projects a complex Python callable onto the scaling basis.
        Usually reached as ScalingProjector(mra, prec, dtype=complex).
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
        "<bra|ket>, conjugating the bra. For the plain integral of the product "
        "use dot(bra.conj(), ket)");
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
        Pointwise complex map through the multiwavelet representation.
        Usually reached as FunctionMap(fmap, prec, dtype=complex).
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

    m.def(
        "build_grid",
        [](FunctionTree<D, ComplexDouble> &out, int scales) { mrcpp::build_grid<D, ComplexDouble>(out, scales); },
        "out"_a,
        "scales"_a);

    m.def(
        "build_grid",
        [](FunctionTree<D, ComplexDouble> &out, std::vector<FunctionTree<D, ComplexDouble> *> &inp, int max_iter) {
            FunctionTreeVector<D, ComplexDouble> vec;
            for (auto *tree : inp) vec.push_back({ComplexDouble(1.0, 0.0), tree});
            mrcpp::build_grid<D, ComplexDouble>(out, vec, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "build_grid",
        [](FunctionTree<D, ComplexDouble> &out,
           std::vector<std::tuple<ComplexDouble, FunctionTree<D, ComplexDouble> *>> &inp,
           int max_iter) {
            FunctionTreeVector<D, ComplexDouble> vec;
            for (auto &t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
            mrcpp::build_grid<D, ComplexDouble>(out, vec, max_iter);
        },
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1);

    m.def(
        "clear_grid",
        [](FunctionTree<D, ComplexDouble> &out) { mrcpp::clear_grid<D, ComplexDouble>(out); },
        "out"_a);

    m.def(
        "refine_grid",
        [](FunctionTree<D, ComplexDouble> &out, int scales) {
            return mrcpp::refine_grid<D, ComplexDouble>(out, scales);
        },
        "out"_a,
        "scales"_a);

    m.def(
        "refine_grid",
        [](FunctionTree<D, ComplexDouble> &out, double prec, bool abs_prec) {
            return mrcpp::refine_grid<D, ComplexDouble>(out, prec, abs_prec);
        },
        "out"_a,
        "prec"_a,
        "abs_prec"_a = false);

    m.def(
        "refine_grid",
        [](FunctionTree<D, ComplexDouble> &out, FunctionTree<D, ComplexDouble> &inp) {
            return mrcpp::refine_grid<D, ComplexDouble>(out, inp);
        },
        "out"_a,
        "inp"_a);

    m.def(
        "power",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           FunctionTree<D, ComplexDouble> &inp,
           double pow,
           int max_iter,
           bool abs_prec) { mrcpp::power<D, ComplexDouble>(prec, out, inp, pow, max_iter, abs_prec); },
        "prec"_a = -1.0,
        "out"_a,
        "inp"_a,
        "pow"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false,
        "Real exponent, principal branch");

    m.def(
        "square",
        [](double prec,
           FunctionTree<D, ComplexDouble> &out,
           FunctionTree<D, ComplexDouble> &inp,
           int max_iter,
           bool abs_prec,
           bool conjugate) { mrcpp::square<D, ComplexDouble>(prec, out, inp, max_iter, abs_prec, conjugate); },
        "prec"_a = -1.0,
        "out"_a,
        "inp"_a,
        "max_iter"_a = -1,
        "abs_prec"_a = false,
        "conjugate"_a = false,
        "conjugate=True gives f * conj(f). The result is real but still a "
        "complex tree, so take .real() to drop the round-off");

    m.def(
        "divergence",
        [](DerivativeOperator<D> &oper, std::vector<FunctionTree<D, ComplexDouble> *> &inp) {
            std::unique_ptr<FunctionTree<D, ComplexDouble>> out{nullptr};
            if (inp.size() == (size_t)D) {
                out = std::make_unique<FunctionTree<D, ComplexDouble>>(inp[0]->getMRA());
                mrcpp::divergence<D, ComplexDouble>(*out, oper, inp);
            }
            return out;
        },
        "oper"_a,
        "inp"_a,
        "No complex gradient in MRCPP; build the components with apply()");
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
}

} // namespace vampyr
