#pragma once

// promote() lives here rather than in complex/complex_bindings.h because the
// operator bindings need it too: a complex kernel cannot be applied to a real
// tree, so the Python __call__ promotes first.

#include <memory>

#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/MWNode.h> // ComplexDouble

namespace vampyr {

/** Deep copy of a real tree into a complex tree (exact, same grid). */
template <int D>
auto promote(mrcpp::FunctionTree<D, double> &inp) -> std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>> {
    return std::unique_ptr<mrcpp::FunctionTree<D, ComplexDouble>>(inp.CopyTreeToComplex());
}

} // namespace vampyr
