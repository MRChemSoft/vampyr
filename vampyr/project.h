/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include <functional>
#include <array>

#include "mrcpp_declarations.h"

namespace vampyr {
template<int D>
void project(double prec, mrcpp::FunctionTree<D> &out,
               std::function<double (std::array<double, D>)> func, int maxIter = -1);
} // namespace vampyr
