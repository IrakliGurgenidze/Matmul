#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include "MatrixUtils.h"

TEST_CASE("Basic Matrix Ground Truth", "[groundTruth]") {
    SECTION("Simple 1:1 match") {
        std::vector<Coord> r1 = {{0, 1}}; // (a=0, b=1)
        std::vector<Coord> r2 = {{1, 2}}; // (b=1, c=2)
        REQUIRE(groundTruthCalc(r1, r2) == 1); // (0,2)
    }

    SECTION("No join") {
        std::vector<Coord> r1 = {{0, 1}, {2, 3}};
        std::vector<Coord> r2 = {{4, 5}, {6, 7}};
        REQUIRE(groundTruthCalc(r1, r2) == 0);
    }

    SECTION("Multiple matches for one b") {
        std::vector<Coord> r1 = {{0, 1}, {2, 1}}; // two a's for b=1
        std::vector<Coord> r2 = {{1, 3}, {1, 4}}; // two c's for b=1
        REQUIRE(groundTruthCalc(r1, r2) == 4); // (0,3), (0,4), (2,3), (2,4)
    }

    SECTION("Duplicate R1 entries should not double count") {
        std::vector<Coord> r1 = {{0, 1}, {0, 1}};
        std::vector<Coord> r2 = {{1, 2}};
        REQUIRE(groundTruthCalc(r1, r2) == 1); // only (0,2)
    }

    SECTION("Multiple distinct b keys") {
        std::vector<Coord> r1 = {{0, 1}, {2, 3}};
        std::vector<Coord> r2 = {{1, 5}, {3, 6}};
        REQUIRE(groundTruthCalc(r1, r2) == 2); // (0,5) and (2,6)
    }
}