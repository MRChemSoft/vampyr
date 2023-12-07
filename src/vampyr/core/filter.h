/*
 *  \date Oct 30, 2023
 *  \author Evgueni Dinvay <evgueni.dinvay@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <MRCPP/core/MWFilter.h>

namespace vampyr {

void filter(pybind11::module &m)
{
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;


    py::class_<MWFilter>(m, "Filter")
        .def(py::init<int, int>(),
             "order"_a,
             "type"_a)
        .def(
            "getFilter",
            &MWFilter::getFilter, py::return_value_policy::reference_internal)
        .def(
            "getReconstructionSubFilter",
            &MWFilter::getReconstructionSubFilter,
            py::return_value_policy::reference_internal,
            py::arg("i"));
}

} // namespace vampyr
