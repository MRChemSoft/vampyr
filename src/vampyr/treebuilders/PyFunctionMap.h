#pragma once

#include <functional>

#include <MRCPP/treebuilders/map.h>

namespace mrcpp {

template <int D> class PyFunctionMap final {
public:
    explicit PyFunctionMap(std::function<double(double)> fmap, double prec = -1.0)
            : precision(prec)
            , func_map(fmap) {}

    std::unique_ptr<FunctionTree<D>> operator()(FunctionTree<D> &inp) {
        // Negative precision will copy grid from input
        auto out = std::make_unique<FunctionTree<D>>(inp.getMRA());
        if (this->precision < 0.0) copy_grid<D>(*out, inp);
        map<D>(this->precision, *out, inp, this->func_map);
        return out;
    }

private:
    double precision;
    std::function<double(double)> func_map;
};

} // namespace mrcpp
