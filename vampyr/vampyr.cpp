/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */


#include "pybind11/pybind11.h"

namespace py = pybind11;

template <int D> void fundamentalTypes(py::module &);
void bases(py::module &);
template<int D> void methods(py::module &);
void pyProject1D(py::module &);
void pyProject2D(py::module &);
void pyProject3D(py::module &);

template<int D> void operators(py::module &);

PYBIND11_MODULE(vampyr, m) {

    pyProject1D(m);
    pyProject2D(m);
    pyProject3D(m);

    fundamentalTypes<1>(m);
    fundamentalTypes<2>(m);
    fundamentalTypes<3>(m);

    bases(m);

    methods<1>(m);
    methods<2>(m);
    methods<3>(m);
    operators<1>(m);
    operators<2>(m);
    operators<3>(m);

}
