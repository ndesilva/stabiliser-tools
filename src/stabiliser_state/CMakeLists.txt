add_library( stabiliser_state 
    stabiliser_state.cpp
    stabiliser_state_from_statevector.cpp
)

target_include_directories( stabiliser_state PUBLIC
    "${PROJECT_SOURCE_DIR}/src/util"
)
