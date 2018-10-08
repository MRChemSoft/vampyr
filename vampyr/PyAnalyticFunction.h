/*
 *  \date May 02, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include "functions/RepresentableFunction.h"
#include <array>

namespace vampyr {

template<int D>
class PyAnalyticFunction : public mrcpp::RepresentableFunction<D> {
public:
    PyAnalyticFunction(std::function<double (std::array<double, D>)> f)
            : func(f) { }

    virtual double evalf(const double *r) const {
        auto val = 0.0;
        auto r_tmp = std::array<double, D>();
        for (auto i = 0; i < D; i++) {
            r_tmp[i] = r[i];
        }
        if (not this->outOfBounds(r)) val = this->func(r_tmp);
        return val;
    }
    virtual double evalf(const std::array<double, D> &r) const {
        return this->evalf(r.data());
    }
protected:
    std::function<double (std::array<double, D>)> func;
};
} // namespace vampyr
