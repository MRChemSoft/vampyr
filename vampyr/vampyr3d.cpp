/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vampyr {
void applys(py::module &);
void bases(py::module &);
void functions(py::module &);
void math(py::module &);
void operators(py::module &);
void operators_grid(py::module &);
void project(py::module &);
void trees(py::module &);
void world(py::module &);

PYBIND11_MODULE(vampyr3d, m) {
    m.doc() = "vampyr3d makes the MRCPP 3D functionality "
              "availible through a python interface";

    applys(m);
    bases(m);
    functions(m);
    math(m);
    operators(m);
    operators_grid(m);
    project(m);
    trees(m);
    world(m);
}
} // namespace vampyr
