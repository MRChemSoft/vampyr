#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"


namespace vampyr {
// Append for FunctionTrees to FunctionTreeVectors
template <int D>
void py_push_back(mrcpp::FunctionTreeVector<D> &vec, float coef, mrcpp::FunctionTree<D> &func_tree) {
    vec.push_back(std::make_tuple(coef, &func_tree));
}

} // namespace vampyr
