/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>

#include <string>

#include "applys/applys.h"

namespace py = pybind11;

namespace vampyr {
void bases(py::module &);
void functions(py::module &);
void math(py::module &);
void operators(py::module &);
void operators_grid(py::module &);
void project(py::module &);
void trees(py::module &);
void world(py::module &);

template <int D> void bind_mr(py::module &mod) noexcept {
    std::string name = "_" + std::to_string(D) + "d";
    py::module sub_mod = mod.def_submodule(name.c_str());

    applys<D>(sub_mod);
    bases(sub_mod);
    functions(sub_mod);
    math(sub_mod);
    operators(sub_mod);
    operators_grid(sub_mod);
    project(sub_mod);
    trees(sub_mod);
    world(sub_mod);
}

PYBIND11_MODULE(_vampyr, m) {
    m.doc() = "VAMPyR makes the MRCPP functionality available through a python interface";

    //bind_mr<1>(m);
    //bind_mr<2>(m);
    bind_mr<3>(m);
}
} // namespace vampyr
