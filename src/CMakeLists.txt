add_subdirectory(f2_helper)
add_subdirectory(stabiliser_state)

set(TEST_LIBRARIES 
    f2_helper_tests
    stabiliser_state_tests
)

add_subdirectory(benchmarking)

add_executable(tests)

target_link_libraries(tests PRIVATE ${TEST_LIBRARIES})
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

include(CTest)
include(Catch)

catch_discover_tests(tests)

add_test( tests COMMAND tests )
