/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <mpi.h>
// see here: https://github.com/mpi4py/mpi4py/issues/19#issuecomment-768143143
#ifdef MSMPI_VER
#define PyMPI_HAVE_MPI_Message 1
#endif
#include <mpi4py/mpi4py.h>
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
#include "treebuilders/maps.h"
#include "treebuilders/project.h"
#include "trees/trees.h"
#include "trees/world.h"

namespace py = pybind11;
using namespace mrcpp;
using namespace pybind11::literals;

struct mpi4py_comm {
    mpi4py_comm() = default;
    mpi4py_comm(MPI_Comm value)
            : value(value) {}
    operator MPI_Comm() { return value; }

    MPI_Comm value;
};

namespace pybind11 {
namespace detail {
template <> struct type_caster<mpi4py_comm> {
public:
    PYBIND11_TYPE_CASTER(mpi4py_comm, _("mpi4py_comm"));

    // Python -> C++
    bool load(handle src, bool) {
        PyObject *py_src = src.ptr();

        // Check that we have been passed an mpi4py communicator
        if (PyObject_TypeCheck(py_src, &PyMPIComm_Type)) {
            // Convert to regular MPI communicator
            value.value = *PyMPIComm_Get(py_src);
        } else {
            return false;
        }

        return !PyErr_Occurred();
    }

    // C++ -> Python
    static handle cast(mpi4py_comm src, return_value_policy /* policy */, handle /* parent */) {
        // Create an mpi4py handle
        return PyMPIComm_New(src.value);
    }
};
} // namespace detail
} // namespace pybind11

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
    advanced_map<D>(sub_mod);
}

template <int D> void bind_vampyr(py::module &mod) noexcept {
    std::string name = "vampyr" + std::to_string(D) + "d";
    py::module sub_mod = mod.def_submodule(name.c_str());

    applys<D>(sub_mod);
    arithmetics<D>(sub_mod);
    project<D>(sub_mod);
    map<D>(sub_mod);

    functions<D>(sub_mod);
    derivatives<D>(sub_mod);
    convolutions<D>(sub_mod);
    trees<D>(sub_mod);
    world<D>(sub_mod);

    bind_advanced<D>(sub_mod);
}

// recieve a communicator and check if it equals MPI_COMM_WORLD
void print_comm(mpi4py_comm comm) {
    if (comm == MPI_COMM_WORLD) {
        std::cout << "Received the world." << std::endl;
    } else {
        std::cout << "Received something else." << std::endl;
    }
}

mpi4py_comm get_comm() {
    return MPI_COMM_WORLD; // Just return MPI_COMM_WORLD for demonstration
}

PYBIND11_MODULE(_vampyr, m) {

    // initialize mpi4py's C-API
    if (import_mpi4py() < 0) {
        // mpi4py calls the Python C API
        // we let pybind11 give us the detailed traceback
        throw py::error_already_set();
    }

    // register the test functions
    m.def("print_comm", &print_comm, "Do something with the mpi4py communicator.");
    m.def("get_comm", &get_comm, "Return some communicator.");

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
