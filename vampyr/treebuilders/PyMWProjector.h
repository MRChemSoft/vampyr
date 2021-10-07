#pragma once

#include <MRCPP/treebuilders/project.h>

namespace mrcpp {

template <int D>
class PyMWProjector final {
public:
    PyMWProjector(const MultiResolutionAnalysis<D> &mra, double prec)
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
