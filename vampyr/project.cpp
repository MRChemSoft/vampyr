/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "treebuilders/project.h"
#include "trees/FunctionTree.h"
#include "PyAnalyticFunction.h"

#include "project.h"
#include <utility>

using namespace mrcpp;

namespace vampyr {
template<int D>
void project(double prec,
               FunctionTree<D> &out,
               std::function<double (std::array<double, D>)> func,
               int maxIter) {
    PyAnalyticFunction<D> inp(std::move(func));
    project(prec, out, inp, maxIter);
}

template<> void project(double prec, FunctionTree<3> &out, std::function<double (std::array<double, 3>)> func, int maxIter);
} // namespace vampyr
