#pragma once

#include <MRCPP/treebuilders/project.h>

namespace mrcpp {

template <int D>
class PyMWProjector final {
public:
    PyMWProjector(const MultiResolutionAnalysis<D> &mra, double prec)
        : min_scale(mra.getRootScale())
        , precision(prec)
        , MRA(mra) {
    }
    PyMWProjector(const MultiResolutionAnalysis<D> &mra, int scale)
        : min_scale(scale)
        , precision(-1.0)
        , MRA(mra) {
    }

    std::unique_ptr<FunctionTree<D>> operator()(RepresentableFunction<D> &func) {
        auto prec = this->precision;
        auto scale = this->min_scale;
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        if (prec > 0.0) {
            // With the adaptive projection we want s+w repr at finest scale
            build_grid<D>(*out, func);
            project<D>(prec, *out, func);
        } else {
            // With the fixed scale projection we want pure s repr at finest scale
            int depth = scale - this->MRA.getRootScale();
            if (depth > 0) {
                // Project at (n-1) and then refine one level to get only s coefs
                build_grid<D>(*out, depth - 1);
                project<D>(-1.0, *out, func);
                refine_grid<D>(*out, 1);
            } else {
                // Need to explicitly clear w coefs after projection to get only s coefs
                project<D>(-1.0, *out, func);
                for (int n = 0; n < out->getNEndNodes(); n++) {
                    int two_d = out->getTDim();
                    int kp1_d = out->getKp1_d();
                    for (int i = 1; i < two_d; i++) {
                        out->getEndMWNode(n).zeroCoefBlock(i, kp1_d);
                    }
                }
            }
        }
        return out;
    }

    std::unique_ptr<FunctionTree<D>> operator()(std::function<double(const Coord<D> &r)> func) {
        auto prec = this->precision;
        auto scale = this->min_scale;
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        if (prec > 0.0) {
            // With the adaptive projection we want s+w repr at finest scale
            project<D>(prec, *out, func);
        } else {
            // With the fixed scale projection we want pure s repr at finest scale
            int depth = scale - this->MRA.getRootScale();
            if (depth > 0) {
                // Project at (n-1) and then refine one level to get only s coefs
                build_grid<D>(*out, depth - 1);
                project<D>(-1.0, *out, func);
                refine_grid<D>(*out, 1);
            } else {
                // Need to explicitly clear w coefs after projection to get only s coefs
                project<D>(-1.0, *out, func);
                for (int n = 0; n < out->getNEndNodes(); n++) {
                    int two_d = out->getTDim();
                    int kp1_d = out->getKp1_d();
                    for (int i = 1; i < two_d; i++) {
                        out->getEndMWNode(n).zeroCoefBlock(i, kp1_d);
                    }
                }
            }
        }
        return out;
    }

private:
    int min_scale;
    double precision;
    MultiResolutionAnalysis<D> MRA;
};

} // namespace mrcpp
