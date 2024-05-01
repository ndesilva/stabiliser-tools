add_library(stabiliser_state_tests OBJECT 
    stabiliser_state.test.cpp 
    stabiliser_state.cpp
    stabiliser_state_from_statevector.test.cpp
    stabiliser_state_from_statevector.cpp
)

target_link_libraries(stabiliser_state_tests PRIVATE Catch2::Catch2WithMain)