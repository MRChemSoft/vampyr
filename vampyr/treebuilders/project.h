#pragma once

#include <pybind11/functional.h>

#include <MRCPP/treebuilders/project.h>

namespace mrcpp {
template <int D>
class MWProjector final {
public:
    MWProjector(const MultiResolutionAnalysis<D> &mra, double prec)
        : precision(prec)
        , MRA(mra) {
    }

    std::unique_ptr<FunctionTree<D>> operator()(RepresentableFunction<D> &func) {
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        build_grid<D>(*out, func);
        project<D>(this->precision, *out, func);
        return out;
    }

    std::unique_ptr<FunctionTree<D>> operator()(std::function<double(const Coord<D> &r)> func) {
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        project<D>(this->precision, *out, func);
        return out;
    }

private:
    double precision;
    MultiResolutionAnalysis<D> MRA;
};
} // namespace mrcpp

namespace vampyr {
template <int D> void project(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<MWProjector<D>>(m, "MWProjector")
        .def(py::init<const MultiResolutionAnalysis<D> &, double>(),
            "mra"_a,
            "prec"_a)
        .def("__call__", [](MWProjector<D> &proj, RepresentableFunction<D> &func){
                return proj(func);
            },
            "func"_a)
        .def("__call__", [](MWProjector<D> &proj, std::function<double (const Coord<D> &r)> func){
                return proj(func);
            },
            "func"_a);

    m.def("project",
          py::overload_cast<double,
                            FunctionTree<D> &,
                            RepresentableFunction<D> &,
                            int,
                            bool>
                            (&mrcpp::project<D>),
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);

    m.def("project", [] (double prec,
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
