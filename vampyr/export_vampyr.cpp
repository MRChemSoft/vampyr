/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/pybind11.h>

#include <string>

#include <MRCPP/version.h>
#include <MRCPP/constants.h>

#include "core/bases.h"
#include "functions/functions.h"
#include "operators/convolutions.h"
#include "operators/derivatives.h"
#include "treebuilders/arithmetics.h"
#include "treebuilders/applys.h"
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

template <int D> void bind_mr(py::module &mod) noexcept {
    applys<D>(mod);
    arithmetics<D>(mod);
    project<D>(mod);
    grids<D>(mod);

    std::string name = "D" + std::to_string(D);
    py::module sub_mod = mod.def_submodule(name.c_str());
    sub_mod.doc() = std::to_string(D) + "-dimensional bindings";

    functions<D>(sub_mod);
    derivatives<D>(sub_mod);
    convolutions<D>(sub_mod);
    trees<D>(sub_mod);
    world<D>(sub_mod);
}

PYBIND11_MODULE(_vampyr, m) {
    using namespace mrcpp;
    m.doc() = "VAMPyR makes the MRCPP functionality available through a Python interface";

    m.def("mrcpp_version", &mrcpp::program_version);

    // Dimension-independent bindings go in the main module
    constants(m);
    bases(m);

    // Dimension-dependent bindings go into submodules
    bind_mr<1>(m);
    bind_mr<2>(m);
    bind_mr<3>(m);
}
} // namespace vampyr
