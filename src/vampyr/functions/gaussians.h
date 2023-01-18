#pragma once

#include <MRCPP/functions/GaussExp.h>
#include <MRCPP/functions/GaussFunc.h>
#include <MRCPP/functions/GaussPoly.h>
#include <MRCPP/functions/function_utils.h>

#include "PyGaussian.h"

namespace vampyr {

template <int D> void gaussians(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // Gaussian class
    py::class_<Gaussian<D>, PyGaussian<D>, RepresentableFunction<D>>(m, "Gaussian",
    R"mydelimiter(
        Parent class to the GaussFunc class.

        Note: This is only a template. Do not use it directly.

    )mydelimiter")
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("pow", py::overload_cast<int>(&Gaussian<D>::getPower, py::const_), "dim"_a)
        .def("exp", py::overload_cast<int>(&Gaussian<D>::getExp, py::const_), "dim"_a = 0)
        .def("pos", &Gaussian<D>::getPos)
        .def("coef", &Gaussian<D>::getCoef)
        .def("periodify", &Gaussian<D>::periodify, "period"_a, "std_dev"_a = 4.0)
        .def("overlap", &Gaussian<D>::calcOverlap, "inp"_a)
        .def("__str__", [](const Gaussian<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        });

    // GaussFunc class
    py::class_<GaussFunc<D>, PyGaussian<D, GaussFunc<D>>, Gaussian<D>>(m, "GaussFunc",
    R"mydelimiter(
        An analytic Gaussian function in d dimensions:

        .. math::

           g(x) = \prod_i^d (x[i] - x_0[i])^{p[i]} \alpha \exp(- \beta \sum_i^d (x[i] - x_0[i])^2 )

        This is a alternative to using a general function definition to project onto the MRA.

        Parameters
        ----------
        beta : scalar
            Exponent of the Gaussian. :math:`\beta`
        alpha : scalar, optional
            Coefficient of Gaussian. :math:`\alpha`
        position : D-dimensional array_like, optional
            Coordinate of the center of Gaussian. :math:`[x_0[0], ..., x_0[d]]`
        poly_exponent: D-dimensional array_like, optional
            Power of the polynomial in front of the Gaussian. :math:`[p[0], ..., p[d]]`

        Returns
        -------
        gaussian : GaussFunc object

    )mydelimiter")
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>(),
             "beta"_a,
             "alpha"_a = 1.0,
             "position"_a = Coord<D>{},
             "poly_exponent"_a = std::array<int, D>{})
        .def(
            "differentiate",
            [](const GaussFunc<D> &gauss, int dir) { return gauss.differentiate(dir).asGaussExp(); },
            "dir"_a, "Differentiates the Gaussian along the specified axis.")
        .def("squaredNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy,
             R"mydelimiter(
             Calculate energy interaction between this Gaussian and an input Gaussian.
             Warning: power has to be a zero vector)mydelimiter");

    // GaussExp class
    py::class_<GaussExp<D>, RepresentableFunction<D>>(m, "GaussExp",
        R"mydelimiter(
       test
        )mydelimiter")
        .def(py::init())
        .def("size",
            py::overload_cast<>(&GaussExp<D>::size, py::const_),
            "Number of Gaussians in the GaussExp")
        .def("func",
            py::overload_cast<int>(&GaussExp<D>::getFunc),
            "term"_a,
            py::return_value_policy::reference_internal)
        .def("append",
            py::overload_cast<const Gaussian<D> &>(&GaussExp<D>::append),
            "Append Gaussians to the end of the GaussExp")
        .def("periodify", &GaussExp<D>::periodify,
            "period"_a,
            "std_dev"_a = 4.0,
            "Make copies of the Gaussian to simulate periodicity, then append it to the GaussExp")
        .def("differentiate",
            &GaussExp<D>::differentiate,
            "dir"_a, "Differentiate all Gaussians in GaussExp along the specified axis")
        .def("squaredNorm", &GaussExp<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussExp<D>::calcCoulombEnergy)
        .def("__str__", [](const GaussExp<D> &func) {
            std::ostringstream os;
            os << func;
            return os.str();
        });
}
} // namespace vampyr
