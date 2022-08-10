#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/Polynomial.h>

namespace vampyr {

class PyScalingFunction : public mrcpp::Polynomial {
public:
    PyScalingFunction(const mrcpp::ScalingBasis &basis, int i, int l, int n)
            : Polynomial(basis.getFunc(i)) {
        this->translate(l);
        this->dilate(std::pow(2, n));
        *this *= std::pow(2, n / 2.0);
    }
};

} // namespace vampyr
