set( SOURCE_FILES
    pauli/pauli.cpp
    stabiliser_state/check_matrix.cpp
    stabiliser_state/stabiliser_state_from_statevector.cpp
    stabiliser_state/stabiliser_state.cpp
)

add_library(fast_stabiliser SHARED ${SOURCE_FILES})
add_library(fast_stabiliser_for_tests ${SOURCE_FILES})

target_include_directories( fast_stabiliser PRIVATE
    "${PROJECT_SOURCE_DIR}/cpp/src"
)

target_include_directories( fast_stabiliser_for_tests PRIVATE
    "${PROJECT_SOURCE_DIR}/cpp/src"
)

set ( PYBIND_FILES
    pauli/pauli_pybind.cpp
)

pybind11_add_module(fast ${PYBIND_FILES})
target_link_libraries(fast PRIVATE fast_stabiliser)