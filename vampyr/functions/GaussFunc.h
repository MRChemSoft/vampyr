#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>

#include <MRCPP/functions/GaussFunc.h>
#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/Gaussian.h>

namespace vampyr {
template <int D> void gauss_func(pybind11::module &m, pybind11::class_<mrcpp::Gaussian<D>> &gaussian) {
    using namespace mrcpp;
    using namespace pybind11::literals;

    pybind11::class_<GaussFunc<D>>(m, "GaussFunc", gaussian)
        .def(pybind11::init<double, double, Coord<D> &, std::array<int, D> &>(),
             "alpha"_a,
             "beta"_a,
             "pos"_a = Coord<D>{},
             "power"_a = std::array<int, D>{})
        .def("evalf", pybind11::overload_cast<const Coord<D> &>(&GaussFunc<D>::evalf, pybind11::const_))
        .def("evalf", pybind11::overload_cast<double, int>(&GaussFunc<D>::evalf, pybind11::const_))
        .def("calcOverlap", pybind11::overload_cast<GaussFunc<D> &>(&GaussFunc<D>::calcOverlap))
        .def("differentiate", &GaussFunc<D>::differentiate)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy)
        .def("getPos", &GaussFunc<D>::getPos)
        .def("getExp", &GaussFunc<D>::getExp)
        .def("getCoef", &GaussFunc<D>::getCoef)
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);
}
} // namespace vampyr
