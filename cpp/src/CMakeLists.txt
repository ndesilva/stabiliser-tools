set( SOURCE_FILES
    pauli/pauli.cpp
    stabiliser_state/check_matrix.cpp
    stabiliser_state/stabiliser_state_from_statevector.cpp
    stabiliser_state/stabiliser_state.cpp
    clifford/clifford.cpp
    clifford/clifford_from_matrix.cpp
)

add_library(fast_stabiliser SHARED ${SOURCE_FILES})
# add_library(fast_stabiliser_for_tests ${SOURCE_FILES})

target_include_directories( fast_stabiliser PRIVATE
    "${PROJECT_SOURCE_DIR}/cpp/src"
)

# target_include_directories( fast_stabiliser_for_tests PRIVATE
#     "${PROJECT_SOURCE_DIR}/cpp/src"
# )

pybind11_add_module(_stab_tools pybind.cpp)
target_include_directories( _stab_tools PRIVATE
    "${PROJECT_SOURCE_DIR}/cpp/src"
)
target_link_libraries(_stab_tools PRIVATE fast_stabiliser)