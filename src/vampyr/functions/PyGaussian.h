#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/Gaussian.h>

#include "PyRepresentableFunction.h"

namespace vampyr {

template <int D, class GaussianBase = mrcpp::Gaussian<D>>
class PyGaussian : public PyRepresentableFunction<D, GaussianBase> {
public:
    /* Inherit the constructors */
    using PyRepresentableFunction<D, GaussianBase>::PyRepresentableFunction;

    GaussianBase *copy() const override { PYBIND11_OVERLOAD_PURE(GaussianBase *, GaussianBase, copy, ); }
    double evalf(const mrcpp::Coord<D> &r) const override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, evalf, r); }
    double evalf1D(double r, int dim) const override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, evalf1D, r, dim); }
    double calcSquareNorm() const override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, calcSquareNorm, ); }
    mrcpp::GaussExp<D> asGaussExp() const override { PYBIND11_OVERLOAD_PURE(mrcpp::GaussExp<D>, GaussianBase, differentiate); }
    mrcpp::GaussPoly<D> differentiate(int dir) const override { PYBIND11_OVERLOAD_PURE(mrcpp::GaussPoly<D>, GaussianBase, differentiate, dir); }
    void setPow(const std::array<int, D> &p) override { PYBIND11_OVERLOAD_PURE(void, GaussianBase, setPow, p); }
    void setPow(int d, int p) override { PYBIND11_OVERLOAD_PURE(void, GaussianBase, setPow, d, p); }
    std::ostream &print(std::ostream &out) const override { PYBIND11_OVERLOAD_PURE(std::ostream &, GaussianBase, print, out); }
};

} // namespace vampyr
