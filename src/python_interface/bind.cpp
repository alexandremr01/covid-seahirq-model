#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "calculate_r0.h"
#include "run_model.h"
#include "parameters.h"

namespace py = pybind11;

PYBIND11_MODULE(cmodels, m) {
    m.doc() = "pybind11 cmodels plugin"; // optional module docstring

    py::class_<Output>(m, "Output")
    .def(py::init<>())
    .def_readwrite("Y_sum", &Output::Y_sum)
    .def_readwrite("YOUT", &Output::YOUT);

    py::class_<PhaseParameters>(m, "PhaseParameters")
    .def(py::init<>())
    .def_readwrite("xI", &PhaseParameters::xI)
    .def_readwrite("beta", &PhaseParameters::beta)
    .def_readwrite("xA", &PhaseParameters::xA);

    py::class_<DynamicParameters>(m, "DynamicParameters")
    .def(py::init<Eigen::Matrix<double, NEA, 1> ,
            Eigen::Matrix<double, NEA, 1> ,
            Eigen::Matrix<double, NEA, NEA> >())
    .def("addPhase", &DynamicParameters::addPhase);

    py::class_<StaticParameters>(m, "StaticParameters")
    .def(py::init<>())
    .def_readwrite("Lambda", &StaticParameters::Lambda)
    .def_readwrite("mu_eq", &StaticParameters::mu_eq)
    .def_readwrite("alpha", &StaticParameters::alpha)
    .def_readwrite("ksi", &StaticParameters::ksi)
    .def_readwrite("rho", &StaticParameters::rho)
    .def_readwrite("phi", &StaticParameters::phi)
    .def_readwrite("eta", &StaticParameters::eta)
    .def_readwrite("a", &StaticParameters::a)
    .def_readwrite("mu_cov", &StaticParameters::mu_cov)
    .def_readwrite("theta", &StaticParameters::theta)
    .def_readwrite("gama", &StaticParameters::gama)
    .def_readwrite("gama_A", &StaticParameters::gama_A)
    .def_readwrite("gama_H", &StaticParameters::gama_H)
    .def_readwrite("gama_HR", &StaticParameters::gama_HR)
    .def_readwrite("gama_RI", &StaticParameters::gama_RI)
    .def_readwrite("gama_RA", &StaticParameters::gama_RA)
    .def_readwrite("gama_RQI", &StaticParameters::gama_RQI)
    .def_readwrite("gama_RQA", &StaticParameters::gama_RQA)
    .def_readwrite("gama_HQI", &StaticParameters::gama_HQI)
    .def_readwrite("Tc", &StaticParameters::Tc)
    .def_readwrite("Tlc", &StaticParameters::Tlc)
    .def_readwrite("rel_pop", &StaticParameters::rel_pop);;

    m.def("model", &run_model, "A function which adds two numbers");
    m.def("calculateR0", &calculateR0, "A function which adds two numbers");

}