add_subdirectory(pauli)
add_subdirectory(stabiliser_state)

target_link_libraries(stabiliser_state PRIVATE pauli)

add_subdirectory(benchmarking)
add_subdirectory(tests)