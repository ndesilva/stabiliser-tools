#include <catch2/catch_test_macros.hpp>

bool placeholder(int number) {
    return 1;
}

TEST_CASE("Empty test case", "[f2 helper]") {
    REQUIRE(placeholder(3) == 1);
}