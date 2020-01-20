/*
 *  \date Mar 16, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */
#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/GaussFunc.h>
#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/Gaussian.h>
#include <MRCPP/functions/RepresentableFunction.h>

template <int D, typename FunctionBase = mrcpp::RepresentableFunction<D>>
class PyRepresentableFunction : public FunctionBase {
public:
    /* Inherit the constructors */
    using FunctionBase::FunctionBase;

    double evalf(const mrcpp::Coord<D> &r) const override { PYBIND11_OVERLOAD_PURE(double, FunctionBase, evalf, r); }
};

// NOTE The following is currently not used
template <int D, class GaussianBase = mrcpp::Gaussian<D>>
class PyGaussian : public PyRepresentableFunction<D, GaussianBase> {
public:
    /* Inherit the constructors */
    using PyRepresentableFunction<D, GaussianBase>::PyRepresentableFunction;

    GaussianBase *copy() const override { PYBIND11_OVERLOAD_PURE(GaussianBase *, GaussianBase, copy, ); }
    double evalf(const mrcpp::Coord<D> &r) const override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, evalf, r); }
    double evalf(double r, int dim) const override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, evalf, r, dim); }
    double calcSquareNorm() override { PYBIND11_OVERLOAD_PURE(double, GaussianBase, calcSquareNorm, ); }
    double calcOverlap(mrcpp::GaussFunc<D> &b) override {
        PYBIND11_OVERLOAD_PURE(double, GaussianBase, calcOverlap, b);
    }
    double calcOverlap(mrcpp::GaussPoly<D> &b) override {
        PYBIND11_OVERLOAD_PURE(double, GaussianBase, calcOverlap, b);
    }
    mrcpp::GaussPoly<D> differentiate(int dir) override {
        PYBIND11_OVERLOAD_PURE(mrcpp::GaussPoly<D>, GaussianBase, differentiate, dir);
    }
    void setPower(const std::array<int, D> &power) override {
        PYBIND11_OVERLOAD_PURE(void, GaussianBase, setPower, power);
    }
    void setPower(int d, int power) override { PYBIND11_OVERLOAD_PURE(void, GaussianBase, setPower, d, power); }
    std::ostream &print(std::ostream &out) const override {
        PYBIND11_OVERLOAD_PURE(std::ostream &, GaussianBase, print, out);
    }
};
