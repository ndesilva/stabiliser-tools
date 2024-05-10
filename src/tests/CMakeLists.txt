find_package( Catch2 3 REQUIRED )

set( SOURCE_FILES
  f2_helper_test.cpp
  stabiliser_state_test.cpp
  stabiliser_state_from_statevector_test.cpp
)

add_executable( tests ${SOURCE_FILES} )

target_include_directories(
  tests PUBLIC
  "${PROJECT_SOURCE_DIR}/src/stabiliser_state"
  "${PROJECT_SOURCE_DIR}/src/util"
)

target_link_libraries( tests PRIVATE Catch2::Catch2WithMain stabiliser_state )

include(CTest)
include(Catch)
catch_discover_tests(tests)