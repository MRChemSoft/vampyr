/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace vampyr {

void operators(py::module &);
void bases(py::module &);
void project(py::module &);
void math(py::module &);
void world(py::module &);
void trees(py::module &);
void applys(py::module &);
void functions(py::module &);

PYBIND11_MODULE(vampyr3d, m) {

    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    bases(m);
    project(m);
    math(m);
    world(m);
    trees(m);
    applys(m);
    operators(m);
    functions(m);

   
}
} // namespace vampyr
