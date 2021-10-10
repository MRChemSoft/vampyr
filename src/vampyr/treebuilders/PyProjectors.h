#pragma once

#include <MRCPP/Printer>
#include <MRCPP/treebuilders/project.h>

namespace mrcpp {

template <int D>
class PyScalingProjector final {
public:
    PyScalingProjector(const MultiResolutionAnalysis<D> &mra, double prec)
        : min_scale(mra.getRootScale())
        , precision(prec)
        , MRA(mra) {
    }
    PyScalingProjector(const MultiResolutionAnalysis<D> &mra, int scale)
        : min_scale(scale)
        , precision(-1.0)
        , MRA(mra) {
            if (this->min_scale < this->MRA.getRootScale()) MSG_ERROR("Invalid scale");
    }

    std::unique_ptr<FunctionTree<D>> operator()(RepresentableFunction<D> &func) {
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        if (this->precision > 0.0) {
            // With the adaptive projection we want s+w repr at finest scale
            build_grid<D>(*out, func);
            project<D>(this->precision, *out, func);
        } else {
            // With the fixed scale projection we want pure s repr at finest scale
            int depth = this->min_scale - this->MRA.getRootScale();
            build_grid<D>(*out, depth);
            project<D>(-1.0, *out, func);

            // Need to explicitly clear w coefs after projection to get only s coefs
            for (int n = 0; n < out->getNEndNodes(); n++) {
                int two_d = out->getTDim();
                int kp1_d = out->getKp1_d();
                for (int i = 1; i < two_d; i++) {
                    out->getEndMWNode(n).zeroCoefBlock(i, kp1_d);
                }
            }
            // Propagate result to root scale
            out->mwTransform(BottomUp);
        }
        return out;
    }

    std::unique_ptr<FunctionTree<D>> operator()(std::function<double(const Coord<D> &r)> func) {
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);
        if (this->precision > 0.0) {
            // With the adaptive projection we want s+w repr at finest scale
            project<D>(this->precision, *out, func);
        } else {
            // With the fixed scale projection we want pure s repr at finest scale
            int depth = this->min_scale - this->MRA.getRootScale();
            build_grid<D>(*out, depth);
            project<D>(-1.0, *out, func);

            // Need to explicitly clear w coefs after projection to get only s coefs
            for (int n = 0; n < out->getNEndNodes(); n++) {
                int two_d = out->getTDim();
                int kp1_d = out->getKp1_d();
                for (int i = 1; i < two_d; i++) {
                    out->getEndMWNode(n).zeroCoefBlock(i, kp1_d);
                }
            }
            // Propagate result to root scale
            out->mwTransform(BottomUp);
        }
        return out;
    }

private:
    int min_scale;
    double precision;
    MultiResolutionAnalysis<D> MRA;
};

template <int D>
class PyWaveletProjector final {
public:
    PyWaveletProjector(const MultiResolutionAnalysis<D> &mra, int scale)
        : min_scale(scale)
        , MRA(mra) {
            if (this->min_scale < this->MRA.getRootScale()) MSG_ERROR("Invalid scale");
    }

    std::unique_ptr<FunctionTree<D>> operator()(RepresentableFunction<D> &func) {
        // With the fixed scale projection we want pure w repr at finest scale
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);

        // Project uniformly at scale n
        int depth = this->min_scale - this->MRA.getRootScale();
        build_grid<D>(*out, depth);
        project<D>(-1.0, *out, func);

        // Need to explicitly clear s coefs after projection to get only w coefs
        for (int n = 0; n < out->getNEndNodes(); n++) {
            int kp1_d = out->getKp1_d();
            out->getEndMWNode(n).zeroCoefBlock(0, kp1_d);
        }

        // Propagate result to root scale
        out->mwTransform(BottomUp);
        return out;
    }

    std::unique_ptr<FunctionTree<D>> operator()(std::function<double(const Coord<D> &r)> func) {
        // With the fixed scale projection we want pure w repr at finest scale
        auto out = std::make_unique<FunctionTree<D>>(this->MRA);

        // Project uniformly at scale n
        int depth = this->min_scale - this->MRA.getRootScale();
        build_grid<D>(*out, depth);
        project<D>(-1.0, *out, func);

        // Need to explicitly clear s coefs after projection to get only w coefs
        for (int n = 0; n < out->getNEndNodes(); n++) {
            int kp1_d = out->getKp1_d();
            out->getEndMWNode(n).zeroCoefBlock(0, kp1_d);
        }

        // Propagate result to root scale
        out->mwTransform(BottomUp);
        return out;
    }

private:
    int min_scale;
    MultiResolutionAnalysis<D> MRA;
};
} // namespace mrcpp
