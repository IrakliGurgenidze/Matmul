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


TEST_CASE("Sparse Matrix Generator", "[generateSparseMatrix]") {
    SECTION("10x10 matrix with 10% sparsity") {
        double sparsity = 0.1;
        int rows = 10, cols = 10;
        int expectedNNZ = static_cast<int>(rows * cols * sparsity);

        auto mat = generateSparseMatrix(sparsity, rows, cols, 42);

        REQUIRE(mat.size() == expectedNNZ);

        // Ensure all entries are within bounds
        for (const auto& coord : mat) {
            REQUIRE(coord.row >= 0);
            REQUIRE(coord.row < rows);
            REQUIRE(coord.col >= 0);
            REQUIRE(coord.col < cols);
        }

        // Ensure there are no duplicate entries
        std::unordered_set<std::pair<int, int>, pair_hash> seen;
        for (const auto& coord : mat) {
            auto [it, inserted] = seen.emplace(coord.row, coord.col);
            REQUIRE(inserted);  // must not be already in the set
        }
    }

    SECTION("100x100 matrix with 1% sparsity, deterministic") {
        double sparsity = 0.01;
        int rows = 100, cols = 100;
        int seed = 123;
        auto mat1 = generateSparseMatrix(sparsity, rows, cols, seed);
        auto mat2 = generateSparseMatrix(sparsity, rows, cols, seed);

        REQUIRE(mat1 == mat2); // deterministic output
    }

    SECTION("Matrix with full density (sparsity = 1.0)") {
        int rows = 5, cols = 5;
        auto mat = generateSparseMatrix(1.0, rows, cols, 34);
        REQUIRE(mat.size() == rows * cols);
    }

    SECTION("Invalid sparsity throws exception") {
        REQUIRE_THROWS(generateSparseMatrix(0.0, 10, 10, 5));
        REQUIRE_THROWS(generateSparseMatrix(1.1, 10, 10, 7));
        REQUIRE_THROWS(generateSparseMatrix(-0.5, 10, 10, 9));
    }
}