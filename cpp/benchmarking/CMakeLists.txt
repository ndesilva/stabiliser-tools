add_subdirectory(stim)

add_executable(stabiliser_benchmark stabiliser_from_statevector_benchmark.cpp)

target_link_libraries(stabiliser_benchmark 
    stim 
    fast_stabiliser
)
target_include_directories(
    stabiliser_benchmark PUBLIC
    "${PROJECT_SOURCE_DIR}/cpp/src"
)