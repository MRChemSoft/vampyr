/*
 *  \date Mar 16, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */
#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/RepresentableFunction.h>

namespace vampyr {

template <int D, typename FunctionBase = mrcpp::RepresentableFunction<D>>
class PyRepresentableFunction : public FunctionBase {
public:
    /* Inherit the constructors */
    using FunctionBase::FunctionBase;

    double evalf(const mrcpp::Coord<D> &r) const override { PYBIND11_OVERLOAD_PURE(double, FunctionBase, evalf, r); }
};

} // namespace vampyr
