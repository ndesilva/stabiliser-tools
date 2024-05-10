add_library( pauli 
    pauli.cpp
)

target_include_directories( pauli PUBLIC
    "${PROJECT_SOURCE_DIR}/src/util"
)