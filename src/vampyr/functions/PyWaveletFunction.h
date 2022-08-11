#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/core/FilterCache.h>
#include <MRCPP/core/MWFilter.h>
#include <MRCPP/functions/Polynomial.h>
#include <MRCPP/functions/RepresentableFunction.h>

namespace vampyr {

class PyWaveletFunction : public mrcpp::RepresentableFunction<1> {
public:
    PyWaveletFunction(const mrcpp::ScalingBasis &basis, int i, int l, int n)
            : index(i)
            , translation(l)
            , scale(n) {

        for (auto i = 0; i <= basis.getScalingOrder(); i++) poly.push_back(basis.getFunc(i));

        auto &lfilters = mrcpp::FilterCache<mrcpp::Legendre>::getInstance();
        auto &ifilters = mrcpp::FilterCache<mrcpp::Interpol>::getInstance();
        int k = basis.getScalingOrder();
        int type = basis.getScalingType();
        switch (type) {
            case mrcpp::Legendre:
                this->filter = &lfilters.get(k);
                break;
            case mrcpp::Interpol:
                this->filter = &ifilters.get(k);
                break;
        }
    }

    double evalf(const mrcpp::Coord<1> &r) const override {

        auto g0 = this->filter->getCompressionSubFilter(2);
        auto g1 = this->filter->getCompressionSubFilter(3);

        double val = 0.0;
        for (auto j = 0; j < poly.size(); j++) {
            auto p0 = poly[j];
            auto p1 = poly[j];
            p0.translate(translation * 2);
            p0.dilate(std::pow(2, scale + 1));
            p0 *= std::pow(2, (scale + 1) / 2.0);
            p1.translate(2 * translation + 1);
            p1.dilate(std::pow(2, scale + 1));
            p1 *= std::pow(2, (scale + 1) / 2.0);
            val += g0(j, index) * p0.evalf(r) + g1(j, index) * p1.evalf(r);
        }

        return val;
    }

private:
    int index, translation, scale;
    mrcpp::MWFilter *filter{nullptr};
    std::vector<mrcpp::Polynomial> poly;
};
} // namespace vampyr
