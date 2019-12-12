#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <iostream>
#include <sstream>

#include <MRCPP/trees/BoundingBox.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
// NOTE These ones will be in a version of MRCPP
namespace detail {
template <typename T> auto stream_collection(const T &coll) -> std::string {
    std::ostringstream os;
    bool first = true;
    os << "[";
    for (auto elem : coll) {
        if (!first) os << ", ";
        os << elem;
        first = false;
    }
    os << "]";
    return os.str();
}
} // namespace detail

template <typename T, size_t D> auto operator<<(std::ostream &os, const std::array<T, D> &coll) -> std::ostream & {
    return (os << detail::stream_collection(coll));
}

template <int D> void world(pybind11::module &m) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    pybind11::class_<BoundingBox<D>>(m, "BoundingBox")
        .def(pybind11::init<int, std::array<int, D>, std::array<int, D>>(),
             "Initiates the BoundingBox of the area functions can exist")
        .def(pybind11::init<std::array<double, D>, bool>(),
             "Initiates the BoundingBox of the area functions can exist, for periodic systems")
        .def("getBoxLength", &BoundingBox<D>::getBoxLength)
        .def("getBoxLengths", &BoundingBox<D>::getBoxLengths)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBounds)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBound)
        .def("getLowerBounds", &BoundingBox<D>::getLowerBounds)
        .def("getLowerBound", &BoundingBox<D>::getLowerBound)
        .def("getUnitLength", &BoundingBox<D>::getUnitLength)
        .def("getScalingFactor", pybind11::overload_cast<int>(&BoundingBox<D>::getScalingFactor, pybind11::const_))
        .def("getCornerIndex", &BoundingBox<D>::getCornerIndex)
        .def("size", pybind11::overload_cast<>(&BoundingBox<D>::size, pybind11::const_))
        .def("size", pybind11::overload_cast<int>(&BoundingBox<D>::size, pybind11::const_))
        .def("getScale", &BoundingBox<D>::getScale);

    pybind11::class_<MultiResolutionAnalysis<D>>(m, "MultiResolutionAnalysis")
        .def(pybind11::init<BoundingBox<D>, ScalingBasis, int>(), "BoundingBox"_a, "Basis Function"_a, "Scale"_a)
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
            os << " unit length       = " << a.getWorldBox().getUnitLengths() << std::endl;
            os << " total boxes       = " << a.getWorldBox().size() << std::endl;
            os << " boxes             = " << boxes << std::endl;
            os << " lower bounds      = " << a.getWorldBox().getLowerBounds() << std::endl;
            os << " upper bounds      = " << a.getWorldBox().getUpperBounds() << std::endl;
            os << " total length      = " << a.getWorldBox().getBoxLengths() << std::endl;
            os << "============================================================" << std::endl;

            return os.str();
        });
}
} // namespace vampyr
