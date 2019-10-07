#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"


#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "trees/MultiResolutionAnalysis.h"
#include "trees/MWNode.h"

#include "utils/mpi_utils.h"

#include "PyRepresentableFunction.h"

#include "operators/DerivativeOperator.h"
#include "operators/GreensKernel.h"

#include "treebuilders/apply.h"
#include "treebuilders/grid.h"

#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

    void abgv_operator(py::module &);
    void bs_operator(py::module &);
    void convolution_operator(py::module &);
    void derivative_kernel(py::module &);
    void derivative_operator(py::module &);
    void greens_kernel(py::module &);
    void helmholtz_kernel(py::module &, py::class_<GreensKernel> &greenskernel);
    void helmholtz_operator(py::module &);
    void identityconvolution_operator(py::module &);
    void identity_kernel(py::module &, py::class_<GreensKernel> &greenskernel);
    void mw_operator(py::module &);
    void operator_state(py::module &);
    void operator_statistics(py::module &);
    void ph_operator(py::module &);
    void poisson_operator(py::module &);


void operators(py::module &m) {

    const auto D = 3;
    
    //GreensKernel
    py::class_<GreensKernel> greenskernel(m, "GreensKernel");

 	  m.def(
        "divergence",
        py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTreeVector<D> &>(&mrcpp::divergence<D>));

    m.def("gradient", py::overload_cast<DerivativeOperator<D> &, FunctionTree<D> &>(&gradient<D>));
 

    abgv_operator(m);
    bs_operator(m);
    convolution_operator(m);
    derivative_kernel(m);
    derivative_operator(m);
    greens_kernel(m);
    helmholtz_kernel(m, greenskernel);
    helmholtz_operator(m);
    identityconvolution_operator(m);
    identity_kernel(m, greenskernel);
    mw_operator(m);
    operator_state(m);
    operator_statistics(m);
    ph_operator(m);
    poisson_operator(m);


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
