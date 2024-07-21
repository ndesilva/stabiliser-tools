find_package( Catch2 3 REQUIRED )

set( TEST_SOURCE_FILES
  f2_helper_test.cpp
  stabiliser_state_test.cpp
  stabiliser_state_from_statevector_test.cpp
  pauli_test.cpp
  check_matrix_test.cpp
  test_util.cpp
  clifford_test.cpp
  clifford_from_matrix_test.cpp
)

add_executable( tests ${TEST_SOURCE_FILES} )

target_include_directories( tests PRIVATE
  "${PROJECT_SOURCE_DIR}/cpp/src"
)

target_link_libraries(tests PRIVATE 
  Catch2::Catch2WithMain 
  fast_stabiliser_for_tests
)

include(CTest)
include(Catch)
catch_discover_tests(tests)