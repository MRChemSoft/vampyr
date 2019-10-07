#include <array>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/MultiResolutionAnalysis.h"


using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void world(py::module &m) {
    const auto D = 3;
    py::class_<BoundingBox<D>>(m, "BoundingBox")
        .def(py::init<int, std::array<int, D>, std::array<int, D>>(),
             "Initiates the BoundingBox of the area functions can exist")
        .def(py::init<std::array<double, D>, bool>(),
             "Initiates the BoundingBox of the area functions can exist, for periodic systems")
        .def("getBoxLength", &BoundingBox<D>::getBoxLength)
        .def("getBoxLengths", &BoundingBox<D>::getBoxLengths)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBounds)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBound)
        .def("getLowerBounds", &BoundingBox<D>::getLowerBounds)
        .def("getLowerBound", &BoundingBox<D>::getLowerBound)
        .def("getUnitLength", &BoundingBox<D>::getUnitLength)
        .def("getScalingFactor", py::overload_cast<int>(&BoundingBox<D>::getScalingFactor, py::const_))
        .def("getCornerIndex", &BoundingBox<D>::getCornerIndex)
        .def("size", py::overload_cast<>(&BoundingBox<D>::size, py::const_))
        .def("size", py::overload_cast<int>(&BoundingBox<D>::size, py::const_))
        .def("getScale", &BoundingBox<D>::getScale);

    py::class_<MultiResolutionAnalysis<D>>(m, "MultiResolutionAnalysis")
        .def(py::init<BoundingBox<D>, ScalingBasis, int>(), "BoundingBox"_a, "Basis Function"_a, "Scale"_a)
        .def("getOrder", &MultiResolutionAnalysis<D>::getOrder, "Returns the order of the scaling basis")
        .def("getWorldBox", &MultiResolutionAnalysis<D>::getWorldBox)
        .def("getScalingBasis", &MultiResolutionAnalysis<D>::getScalingBasis)
        .def("getMaxDepth", &MultiResolutionAnalysis<D>::getMaxDepth)
        .def("getMaxScale", &MultiResolutionAnalysis<D>::getMaxScale)
        .def("print", [](const MultiResolutionAnalysis<D> &a) {
            auto boxes = std::array<int, D>{a.getWorldBox().size(0), a.getWorldBox().size(1), a.getWorldBox().size(2)};
            py::print("============================================================");
            py::print("                  MultiResolution Analysis                  ");
            py::print("------------------------------------------------------------");
            py::print(" scaling order     = ", a.getScalingBasis().getScalingOrder());
            if (a.getScalingBasis().getScalingType() == Legendre) {
                py::print(" polynomial type   =  Legendre polynomial");
            } else if (a.getScalingBasis().getScalingType() == Interpol) {
                py::print(" polynomial type   =  Interpolating polynomial");
            } else {
                py::print(" polynomial type   =  Unknown");
            }
            py::print("------------------------------------------------------------");
            py::print(" unit length       = ", a.getWorldBox().getUnitLengths());
            py::print(" total boxes       = ", a.getWorldBox().size());
            py::print(" boxes             = ", boxes);
            py::print(" lower bounds      = ", a.getWorldBox().getLowerBounds());
            py::print(" upper bounds      = ", a.getWorldBox().getUpperBounds());
            py::print(" total length      = ", a.getWorldBox().getBoxLengths());
            py::print("============================================================");
        });

}
} // namespace vampyr
