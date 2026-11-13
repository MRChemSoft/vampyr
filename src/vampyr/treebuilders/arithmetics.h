#pragma once

#include <pybind11/complex.h>

#include <MRCPP/treebuilders/add.h>
#include <MRCPP/treebuilders/multiply.h>

namespace vampyr {
template <int D> void arithmetics(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "sum",
        [](std::vector<FunctionTree<D, double> *> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, double>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = inp[0]->getMRA();
                out = std::make_unique<FunctionTree<D, double>>(mra);
                FunctionTreeVector<D, double> vec;
                for (auto* tree : inp) vec.push_back({1.0, tree});
                build_grid(*out, vec);
                add(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def(
        "sum",
        [](std::vector<std::tuple<double, FunctionTree<D, double> *>> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, double>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = std::get<1>(inp[0])->getMRA();
                out = std::make_unique<FunctionTree<D, double>>(mra);
                FunctionTreeVector<D, double> vec;
                for (auto& t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
                build_grid(*out, vec);
                add(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def("dot", [](FunctionTree<D, double> &bra, FunctionTree<D, double> &ket) { return mrcpp::dot<D, double>(bra, ket); }, "bra"_a, "ket"_a);

    m.def(
        "dot",
        [](std::vector<FunctionTree<D, double> *> &inp_a, std::vector<FunctionTree<D, double> *> &inp_b) {
            auto out = std::unique_ptr<FunctionTree<D, double>>(nullptr);
            if ((inp_a.size() > 0) && (inp_b.size() == inp_a.size())) {
                auto &mra = inp_a[0]->getMRA();
                out = std::make_unique<FunctionTree<D, double>>(mra);
                auto out_vec = FunctionTreeVector<D, double>();
                for (size_t i = 0; i < inp_a.size(); ++i) {
                    auto *prod = new FunctionTree<D, double>(mra);
                    build_grid(*prod, *inp_a[i]);
                    build_grid(*prod, *inp_b[i]);
                    build_grid(*prod, 1);
                    multiply(-1.0, *prod, 1.0, *inp_a[i], *inp_b[i]);
                    out_vec.push_back({1.0, prod});
                }
                build_grid(*out, out_vec); // Union grid
                add(-1.0, *out, out_vec);
                clear(out_vec, true);
            }
            return out;
        },
        "inp_a"_a,
        "inp_b"_a);

    m.def(
        "prod",
        [](std::vector<FunctionTree<D, double> *> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, double>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = inp[0]->getMRA();
                out = std::make_unique<FunctionTree<D, double>>(mra);
                FunctionTreeVector<D, double> vec;
                for (auto* tree : inp) vec.push_back({1.0, tree});
                build_grid(*out, vec); // Union grid
                build_grid(*out, 1);   // One extra refinement
                multiply(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def(
        "prod",
        [](std::vector<std::tuple<double, FunctionTree<D, double> *>> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, double>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = std::get<1>(inp[0])->getMRA();
                out = std::make_unique<FunctionTree<D, double>>(mra);
                FunctionTreeVector<D, double> vec;
                for (auto& t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
                build_grid(*out, vec); // Union grid
                build_grid(*out, 1);   // One extra refinement
                multiply(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

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

    m.def(
        "sum",
        [](std::vector<FunctionTree<D, ComplexDouble> *> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, ComplexDouble>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = inp[0]->getMRA();
                out = std::make_unique<FunctionTree<D, ComplexDouble>>(mra);
                FunctionTreeVector<D, ComplexDouble> vec;
                for (auto *tree : inp) vec.push_back({1.0, tree});
                build_grid(*out, vec);
                add(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def(
        "sum",
        [](std::vector<std::tuple<ComplexDouble, FunctionTree<D, ComplexDouble> *>> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, ComplexDouble>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = std::get<1>(inp[0])->getMRA();
                out = std::make_unique<FunctionTree<D, ComplexDouble>>(mra);
                FunctionTreeVector<D, ComplexDouble> vec;
                for (auto &t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
                build_grid(*out, vec);
                add(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def(
        "prod",
        [](std::vector<FunctionTree<D, ComplexDouble> *> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, ComplexDouble>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = inp[0]->getMRA();
                out = std::make_unique<FunctionTree<D, ComplexDouble>>(mra);
                FunctionTreeVector<D, ComplexDouble> vec;
                for (auto *tree : inp) vec.push_back({1.0, tree});
                build_grid(*out, vec); // Union grid
                build_grid(*out, 1);   // One extra refinement
                multiply(-1.0, *out, vec);
            }
            return out;
        },
        "inp"_a);

    m.def(
        "prod",
        [](std::vector<std::tuple<ComplexDouble, FunctionTree<D, ComplexDouble> *>> &inp) {
            auto out = std::unique_ptr<FunctionTree<D, ComplexDouble>>(nullptr);
            if (inp.size() > 0) {
                auto &mra = std::get<1>(inp[0])->getMRA();
                out = std::make_unique<FunctionTree<D, ComplexDouble>>(mra);
                FunctionTreeVector<D, ComplexDouble> vec;
                for (auto &t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
                build_grid(*out, vec); // Union grid
                build_grid(*out, 1);   // One extra refinement
                multiply(-1.0, *out, vec);
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
          [](double prec, FunctionTree<D, double> &out, double a, FunctionTree<D, double> &inp_a, double b, FunctionTree<D, double> &inp_b, int max_iter, bool abs_prec) {
              mrcpp::add<D, double>(prec, out, a, inp_a, b, inp_b, max_iter, abs_prec);
          },
          "prec"_a = -1.0,
          "out"_a,
          "a"_a = 1.0,
          "inp_a"_a,
          "b"_a = 1.0,
          "inp_b"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("add",
          [](double prec, FunctionTree<D, double> &out, std::vector<FunctionTree<D, double> *> &inp, int max_iter, bool abs_prec) {
              FunctionTreeVector<D, double> vec;
              for (auto* tree : inp) vec.push_back({1.0, tree});
              mrcpp::add<D, double>(prec, out, vec, max_iter, abs_prec);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("add",
          [](double prec, FunctionTree<D, double> &out, std::vector<std::tuple<double, FunctionTree<D, double> *>> &inp, int max_iter, bool abs_prec) {
              FunctionTreeVector<D, double> vec;
              for (auto& t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
              mrcpp::add<D, double>(prec, out, vec, max_iter, abs_prec);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("multiply",
          [](double prec, FunctionTree<D, double> &out, double c, FunctionTree<D, double> &inp_a, FunctionTree<D, double> &inp_b, int max_iter, bool abs_prec, bool use_max_norms) {
              mrcpp::multiply<D, double>(prec, out, c, inp_a, inp_b, max_iter, abs_prec, use_max_norms);
          },
          "prec"_a = -1.0,
          "out"_a,
          "c"_a = 1.0,
          "inp_a"_a,
          "inp_b"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def("multiply",
          [](double prec, FunctionTree<D, double> &out, std::vector<FunctionTree<D, double> *> &inp, int max_iter, bool abs_prec, bool use_max_norms) {
              FunctionTreeVector<D, double> vec;
              for (auto* tree : inp) vec.push_back({1.0, tree});
              mrcpp::multiply<D, double>(prec, out, vec, max_iter, abs_prec, use_max_norms);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def("multiply",
          [](double prec, FunctionTree<D, double> &out, std::vector<std::tuple<double, FunctionTree<D, double> *>> &inp, int max_iter, bool abs_prec, bool use_max_norms) {
              FunctionTreeVector<D, double> vec;
              for (auto& t : inp) vec.push_back({std::get<0>(t), std::get<1>(t)});
              mrcpp::multiply<D, double>(prec, out, vec, max_iter, abs_prec, use_max_norms);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false,
          "use_max_norms"_a = false);

    m.def(
        "dot",
        [](double prec, FunctionTree<D, double> &out, std::vector<FunctionTree<D, double>*> &inp_a, std::vector<FunctionTree<D, double>*> &inp_b, int maxIter, bool abs_prec) {
            FunctionTreeVector<D, double> vec_a, vec_b;
            for (auto* t : inp_a) vec_a.push_back({1.0, t});
            for (auto* t : inp_b) vec_b.push_back({1.0, t});
            mrcpp::dot<D, double>(prec, out, vec_a, vec_b, maxIter, abs_prec);
        },
        "prec"_a,
        "out"_a,
        "inp_a"_a,
        "inp_b"_a,
        "maxIter"_a = -1,
        "abs_prec"_a = false);

    m.def("power",
          [](double prec, FunctionTree<D, double> &out, FunctionTree<D, double> &inp, double pow, int max_iter, bool abs_prec) {
              mrcpp::power<D, double>(prec, out, inp, pow, max_iter, abs_prec);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "pow"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("square",
          [](double prec, FunctionTree<D, double> &out, FunctionTree<D, double> &inp, int max_iter, bool abs_prec) {
              mrcpp::square<D, double>(prec, out, inp, max_iter, abs_prec);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

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
}
} // namespace vampyr
