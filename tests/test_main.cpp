#define CATCH_CONFIG_MAIN // This tells Catch2 to provide a main() function
#include <catch2/catch_test_macros.hpp>
#include <cstdint>

// Example Function
int add(int a, int b) {
    return a + b;
}

TEST_CASE("Addition works correctly", "[math]") {
    REQUIRE(add(1, 1) == 2);
    REQUIRE(add(-1, -1) == -2);
    REQUIRE(add(5, -5) == 0);
}