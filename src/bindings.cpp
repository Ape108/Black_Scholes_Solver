#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // autoconverts std::vector to python lists
#include "black_scholes.hpp"

namespace py = pybind11;

PYBIND11_MODULE(black_scholes_solver, m) {
    m.doc() = "C++ Black-Scholes Finite Difference Solver";

    // Expose the GridParams Struct to python
    py::class_<GridParams>(m, "GridParams")
        .def(py::init<>()) // Allow instantiation in python: grid = black_scholes_solver.GridParams()
        .def_readwrite("price_ceiling", &GridParams::price_ceiling)
        .def_readwrite("time_to_maturity", &GridParams::time_to_maturity)
        .def_readwrite("num_price_steps", &GridParams::num_price_steps)
        .def_readwrite("num_time_steps", &GridParams::num_time_steps);

    // Expose the MarketParams Struct
    py::class_<MarketParams>(m, "MarketParams")
        .def(py::init<>())
        .def_readwrite("volatility", &MarketParams::volatility)
        .def_readwrite("risk_free_interest", &MarketParams::risk_free_interest)
        .def_readwrite("strike_price", &MarketParams::strike_price);

    // Bind core formulation function
    m.def("formulate_black_scholes", &formulate_black_scholes,
          "Solves the Black-Scholes PDE using Crank-Nicolson",
          py::arg("grid"), py::arg("market"));
}