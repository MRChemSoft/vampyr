#pragma once

#include <memory>

#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h> // ComplexDouble

namespace vampyr {

/** Copy of a real tree as a complex tree on the same grid. */
template <int D>
auto promote(mrcpp::FunctionTree<D, double> &inp) -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    return std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>>(inp.CopyTreeToComplex());
}

} // namespace vampyr
