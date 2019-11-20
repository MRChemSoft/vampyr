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
        .def("__repr__", [](const MultiResolutionAnalysis<D> &a) {
            std::ostringstream os;
                
            auto boxes = std::array<int, D>{a.getWorldBox().size(0), a.getWorldBox().size(1), a.getWorldBox().size(2)};
            os << "============================================================" << std::endl; 
            os << "                  MultiResolution Analysis                  " << std::endl; 
            os << "------------------------------------------------------------" << std::endl; 
            os << " scaling order     = " << a.getScalingBasis().getScalingOrder() << std::endl; 
            if (a.getScalingBasis().getScalingType() == Legendre) { 
                os << " polynomial type   =  Legendre polynomial" << std::endl; 
            } else if (a.getScalingBasis().getScalingType() == Interpol) {
                os << " polynomial type   =  Interpolating polynomial" << std::endl; 
            } else {
                os << " polynomial type   =  Unknown" << std::endl; 
            }
            os << "------------------------------------------------------------" << std::endl;
            os << " unit length       = " << a.getWorldBox().getUnitLengths()[0] << ", " << a.getWorldBox().getUnitLengths()[1] << ", " << a.getWorldBox().getUnitLengths()[2] << std::endl;
            os << " total boxes       = " << a.getWorldBox().size() << std::endl;
            os << " boxes             = " << boxes [0] << ", " << boxes [1] << ", " << boxes [2] << std::endl;
            os << " lower bounds      = " << a.getWorldBox().getLowerBounds()[0] << ", " << a.getWorldBox().getLowerBounds()[1] << ", " << a.getWorldBox().getLowerBounds()[2] << std::endl;
            os << " upper bounds      = " << a.getWorldBox().getUpperBounds()[0] << ", " << a.getWorldBox().getUpperBounds()[1] << ", " << a.getWorldBox().getUpperBounds()[2] << std::endl;
            os << " total length      = " << a.getWorldBox().getBoxLengths()[0] << ", " << a.getWorldBox().getBoxLengths()[1] << ", " << a.getWorldBox().getBoxLengths()[2] << std::endl;
            os << "============================================================" << std::endl;
       
            return os.str();


        });

}
} // namespace vampyr
