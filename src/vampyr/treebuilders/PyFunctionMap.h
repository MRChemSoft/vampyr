#pragma once

#include <functional>

#include <MRCPP/treebuilders/treeMap.h>

namespace mrcpp {

template <int D, typename T = double> class PyFunctionMap final {
public:
    explicit PyFunctionMap(std::function<T(T)> fmap, double prec = -1.0)
            : precision(prec)
            , func_map(fmap) {}

    std::unique_ptr<FunctionTree<D, T>> operator()(FunctionTree<D, T> &inp) {
        // Negative precision will copy grid from input
        auto out = std::make_unique<FunctionTree<D, T>>(inp.getMRA());
        if (this->precision < 0.0) copy_grid<D, T>(*out, inp);
        treeMap<D, T>(this->precision, *out, inp, this->func_map);
        return out;
    }

private:
    double precision;
    std::function<T(T)> func_map;
};

} // namespace mrcpp
