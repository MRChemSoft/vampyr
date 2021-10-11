#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <array>
#include <iostream>
#include <sstream>

#include <MRCPP/trees/BoundingBox.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {

template <int D> void world(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    pybind11::class_<BoundingBox<D>>(m, "BoundingBox")
        .def(py::init<int,
                      std::array<int, D> &,
                      std::array<int, D> &,
                      std::array<double, D> &,
                      bool>(),
                      "scale"_a = 0,
                      "corner"_a = std::array<int, D>{},
                      "nboxes"_a = std::array<int, D>{},
                      "scaling"_a = std::array<double, D>{},
                      "pbc"_a = false)
        .def("boxLengths", &BoundingBox<D>::getBoxLengths)
        .def("boxLength", &BoundingBox<D>::getBoxLength, "dim"_a)
        .def("upperBounds", &BoundingBox<D>::getUpperBounds)
        .def("upperBound", &BoundingBox<D>::getUpperBound, "dim"_a)
        .def("lowerBounds", &BoundingBox<D>::getLowerBounds)
        .def("lowerBound", &BoundingBox<D>::getLowerBound, "dim"_a)
        .def("unitLengths", &BoundingBox<D>::getUnitLengths)
        .def("unitLength", &BoundingBox<D>::getUnitLength, "dim"_a)
        .def("scalingFactors", &BoundingBox<D>::getScalingFactors)
        .def("scalingFactor", &BoundingBox<D>::getScalingFactor, "dim"_a)
        .def("isPeriodic", &BoundingBox<D>::isPeriodic)
        .def("scale", &BoundingBox<D>::getScale)
        .def("size", pybind11::overload_cast<>(&BoundingBox<D>::size, pybind11::const_))
        .def("size", pybind11::overload_cast<int>(&BoundingBox<D>::size, pybind11::const_), "dim"_a)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__str__", [](const BoundingBox<D> &box) {
            std::ostringstream os;
            os << box;
            return os.str();
        });

    pybind11::class_<MultiResolutionAnalysis<D>>(m, "MultiResolutionAnalysis")
        .def(py::init<std::array<int, 2>,
                      int,
                      int>(),
                      "box"_a,
                      "order"_a,
                      "max_depth"_a = 30)
        .def(py::init<const BoundingBox<D>,
                      int,
                      int>(),
                      "box"_a,
                      "order"_a,
                      "max_depth"_a = 30)
        .def(py::init<const BoundingBox<D>,
                      const ScalingBasis,
                      int>(),
                      "box"_a,
                      "basis"_a,
                      "max_depth"_a = 30)
        .def("basis", &MultiResolutionAnalysis<D>::getScalingBasis)
        .def("world", &MultiResolutionAnalysis<D>::getWorldBox)
        .def("maxDepth", &MultiResolutionAnalysis<D>::getMaxDepth)
        .def("maxScale", &MultiResolutionAnalysis<D>::getMaxScale)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__str__", [](const MultiResolutionAnalysis<D> &mra) {
            std::ostringstream os;
            os << "================================================================" << std::endl;
            os << "                    MultiResolution Analysis                    " << std::endl;
            os << "----------------------------------------------------------------" << std::endl;
            os << mra.getScalingBasis()                                              << std::endl;
            os << "----------------------------------------------------------------" << std::endl;
            os << mra.getWorldBox()                                                  << std::endl;
            os << "================================================================" << std::endl;
            return os.str();
        });
}
} // namespace vampyr
