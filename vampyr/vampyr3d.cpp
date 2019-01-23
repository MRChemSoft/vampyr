/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "trees/BoundingBox.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "trees/MWTree.h"
#include "trees/MultiResolutionAnalysis.h"

#include "PyRepresentableFunction.h"

#include "operators/ABGVOperator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "operators/HelmholtzOperator.h"
#include "operators/PoissonOperator.h"

#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

#include "functions/GaussFunc.h"
#include "functions/GaussPoly.h"
#include "functions/Gaussian.h"
#include "functions/Polynomial.h"
#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void bases(py::module &);

PYBIND11_MODULE(vampyr3d, m) {

    const auto D = 3;

    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    bases(m);

    m.def(
        "project",
        py::overload_cast<double, FunctionTree<D> &, std::function<double(const Coord<D> &)>, int>(&mrcpp::project<D>),
        "precision"_a,
        "output_tree"_a,
        "function"_a,
        "maxIter"_a = -1,
        "Projects an analytic function onto a FunctionTree");

    m.def("project",
          py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int>(&mrcpp::project<D>),
          "precision"_a,
          "output_tree"_a,
          "GaussFunc"_a,
          "maxIter"_a = -1);

    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, "RepresentableFunction");
    repfunc.def(py::init<>());

    py::class_<RepresentableFunction<1>, PyRepresentableFunction<1>> repfunc1d(m, "RepresentableFunction1D");
    repfunc.def(py::init<>());

    py::class_<Gaussian<D>> gaussian(m, "Gaussian", repfunc);

    py::class_<GaussFunc<D>>(m, "GaussFunc", gaussian)
        .def(py::init<double, double, Coord<D> &, std::array<int, D> &>())
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussFunc<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<double, int>(&GaussFunc<D>::evalf, py::const_))
        .def("calcSquareNorm", &GaussFunc<D>::calcSquareNorm)
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);

    py::class_<GaussPoly<D>>(m, "GaussPoly", gaussian)
        .def(py::init<double, double, const Coord<D> &, const std::array<int, D> &>())
        .def("setPoly", &GaussPoly<D>::setPoly)
        .def("differentiate", &GaussPoly<D>::differentiate)
        .def("evalf", py::overload_cast<const Coord<D> &>(&GaussPoly<D>::evalf, py::const_))
        .def("evalf", py::overload_cast<const double, int>(&GaussPoly<D>::evalf, py::const_));

    py::class_<Polynomial>(m, "Polynomial", repfunc1d)
        .def(py::init<int, const double *, const double *>(), "order"_a, "a"_a = nullptr, "b"_a = nullptr)
        .def("evalf", py::overload_cast<double>(&Polynomial::evalf, py::const_))
        .def("setCoefs", &Polynomial::setCoefs)
        .def("getCoefs", py::overload_cast<>(&Polynomial::getCoefs))
        .def("size", &Polynomial::size);

    py::class_<BoundingBox<D>>(m, "BoundingBox")
        .def(py::init<int, std::array<int, D>, std::array<int, D>>(),
             "Initiates the BoundingBox of the area functions can exist")
        .def("getBoxLength", &BoundingBox<D>::getBoxLength)
        .def("getBoxLengths", &BoundingBox<D>::getBoxLengths)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBounds)
        .def("getUpperBounds", &BoundingBox<D>::getUpperBound)
        .def("getLowerBounds", &BoundingBox<D>::getLowerBounds)
        .def("getLowerBound", &BoundingBox<D>::getLowerBound)
        .def("getUnitLength", &BoundingBox<D>::getUnitLength)
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

    py::class_<MWTree<D>> mwtree(m, "MWTree");
    mwtree.def(py::init<MultiResolutionAnalysis<D>>()).def("getSquareNorm", &MWTree<D>::getSquareNorm);

    py::class_<FunctionTree<D>> functree(m, "FunctionTree", mwtree);
    functree.def(py::init<MultiResolutionAnalysis<D>>())
        .def("integrate", &FunctionTree<D>::integrate, "Integral of the FunctionTree over all space")
        .def("clear", &FunctionTree<D>::clear, "Clears the FunctionTree")
        .def("normalize", &FunctionTree<D>::normalize, "Rescale the function by its norm, fixed grid")
        .def("rescale", &FunctionTree<D>::rescale, "Rescales the function")
        .def("evalf",
             py::overload_cast<const Coord<D> &>(&FunctionTree<D>::evalf),
             "Returns the function value at a given point");

    m.def("apply",
          py::overload_cast<double, FunctionTree<D> &, ConvolutionOperator<D> &, FunctionTree<D> &, int>(&apply<D>),
          "precision"_a,
          "output_tree"_a,
          "oper"_a,
          "input_tree"_a,
          "maxIter"_a = -1,
          "Applies a given convolution operator onto a FunctionTree");

    m.def("apply",
          py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTree<D> &, int>(&apply<D>),
          "out_tree"_a,
          "DerivativeOperator"_a,
          "input_tree"_a,
          "dir"_a = -1,
          "Applies a given derivative operator onto a FunctionTree");

    m.def("add",
          py::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, double, FunctionTree<D> &, int>(
              &add<D>),
          "precision"_a,
          "output_tree"_a,
          "a"_a,
          "tree_a"_a,
          "b"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "Addition: output_tree = a*tree_a + b*tree_b");

    m.def("add",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int>(&add<D>),
          "precision"_a,
          "output_tree"_a,
          "FunctionTreeVector"_a,
          "maxIter"_a = -1,
          "Addition: output_tree = sum_i c_i f_i(x)");

    m.def("multiply",
          py::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, FunctionTree<D> &, int>(&multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "c"_a,
          "tree_a"_a,
          "tree_b"_a,
          "maxIter"_a = -1,
          "Multiplication: output_tree = c*tree_a*tree_b");

    m.def("multiply",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int>(&multiply<D>),
          "precision"_a,
          "output_tree"_a,
          "inpit_vector"_a,
          "maxIter"_a = -1,
          "Multiplication: output_tree = prod_i c_i f_i(x)");

    m.def("dot", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&dot<D>), "bra"_a, "ket"_a);

    m.def("dot",
          py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, FunctionTreeVector<D> &, int>(&dot<D>),
          "prec"_a,
          "out"_a,
          "inp_a"_a,
          "inp_b"_a,
          "maxIter"_a = -1);

    m.def("power",
          py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, double, int>(&power<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "power"_a,
          "maxIter"_a = -1);

    m.def("square",
          py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, int>(&square<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("map",
          py::overload_cast<double, FunctionTree<D> &, FunctionTree<D> &, RepresentableFunction<D> &>(&map<D>),
          "prec"_a,
          "out"_a,
          "inp"_a,
          "Func"_a);

    py::class_<DerivativeOperator<D>> deriv(m, "Derivative Operator");
    deriv.def(py::init<MultiResolutionAnalysis<D>>());

    py::class_<ABGVOperator<D>>(m, "ABGVOperator", deriv)
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "Derivative Operator: The ABGVOperator for differentiation");

    m.def(
        "divergence",
        py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTreeVector<D> &>(&mrcpp::divergence<D>));

    m.def("gradient", py::overload_cast<DerivativeOperator<D> &, FunctionTree<D> &>(&gradient<D>));

    py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
    convop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    py::class_<PoissonOperator>(m, "PoissonOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, double>(),
             "MRA"_a,
             "precision"_a,
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");

    py::class_<HelmholtzOperator>(m, "HelmholtzOperator", convop)
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "MRA"_a,
             "mu"_a,
             "precision"_a,
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, const RepresentableFunction<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTree<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("build_grid",
          py::overload_cast<FunctionTree<D> &, FunctionTreeVector<D> &, int>(&build_grid<D>),
          "out"_a,
          "inp"_a,
          "maxIter"_a = -1);

    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_grid<D>));
    m.def("copy_func", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&copy_func<D>));

    m.def("clear_grid", py::overload_cast<FunctionTree<D> &>(&clear_grid<D>));
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, int>(&refine_grid<D>), "out"_a, "scales"_a);

    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, double>(&refine_grid<D>), "out"_a, "prec"_a);
    m.def("refine_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &>(&refine_grid<D>), "out"_a, "inp"_a);
}
} // namespace vampyr
