/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>

#include <string>

#include <MRCPP/constants.h>
#include <MRCPP/version.h>

#include "core/bases.h"
#include "functions/functions.h"
#include "operators/convolutions.h"
#include "operators/derivatives.h"
#include "treebuilders/applys.h"
#include "treebuilders/arithmetics.h"
#include "treebuilders/grids.h"
#include "treebuilders/project.h"
#include "trees/trees.h"
#include "trees/world.h"

namespace py = pybind11;
using namespace mrcpp;
using namespace pybind11::literals;

namespace vampyr {

void constants(py::module &m) {
    py::enum_<Traverse>(m, "Traverse")
        .value("TopDown", Traverse::TopDown)
        .value("BottomUp", Traverse::BottomUp)
        .export_values();
    py::enum_<Iterator>(m, "Iterator")
        .value("Lebesgue", Iterator::Lebesgue)
        .value("Hilbert", Iterator::Hilbert)
        .export_values();
}

template <int D> void bind_advanced(py::module &mod) noexcept {
    py::module sub_mod = mod.def_submodule("advanced");

    advanced_applys<D>(sub_mod);
    advanced_arithmetics<D>(sub_mod);
    advanced_project<D>(sub_mod);
    advanced_grids<D>(sub_mod);
}

template <int D> void bind_vampyr(py::module &mod) noexcept {
    std::string name = "vampyr" + std::to_string(D) + "d";
    py::module sub_mod = mod.def_submodule(name.c_str());

    applys<D>(sub_mod);
    arithmetics<D>(sub_mod);
    project<D>(sub_mod);

    functions<D>(sub_mod);
    derivatives<D>(sub_mod);
    convolutions<D>(sub_mod);
    trees<D>(sub_mod);
    world<D>(sub_mod);

    bind_advanced<D>(sub_mod);
}

PYBIND11_MODULE(_vampyr, m) {
    m.doc() = R"pbdoc(
        VAMPyR
        ------

        VAMPyR makes the MRCPP functionality available through a Python interface.

        .. currentmodule:: vampyr

        .. autosummary::
           :toctree: generate
    )pbdoc";

    m.attr("__version__") = VERSION_INFO;

    m.def("mrcpp_version", &program_version, "Return version of the underlying MRCPP library.");

    // Dimension-independent bindings go in the main module
    constants(m);

    // Dimension-dependent bindings go into submodules
    bind_vampyr<1>(m);
    bind_vampyr<2>(m);
    bind_vampyr<3>(m);
    bases(m);
}
} // namespace vampyr
