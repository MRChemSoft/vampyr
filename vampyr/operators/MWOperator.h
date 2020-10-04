#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/MWOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
void mw_operator(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    //MWOperator
    py::class_<MWOperator>(m, "MWOperator")
        .def(py::init<const MultiResolutionAnalysis<2> &>())
        .def("size", &MWOperator::size, "Size")
        .def("push_back", &MWOperator::push_back, "Push")
        .def("clear", &MWOperator::clear, "Clear")
        .def("getMaxBandWidth", &MWOperator::getMaxBandWidth, "Get")
        .def("calcBandWidths", &MWOperator::calcBandWidths, "Calc")
        .def("clearBandWidths", &MWOperator::clearBandWidths, "Clear");
}
} // namespace vampyr
