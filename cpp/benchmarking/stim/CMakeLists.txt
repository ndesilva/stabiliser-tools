add_library(stim
    circuit_vs_amplitudes.cc
    twiddle.h
    vector_simulator.cc
)

pybind11_add_module(stim_mock stim_pybind.cpp)

target_link_libraries(stim_mock PRIVATE stim)