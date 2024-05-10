include_directories("${PROJECT_SOURCE_DIR}/src/f2_helper")

add_library( stabiliser_state 
    stabiliser_state.cpp
    stabiliser_state_from_statevector.cpp
)
