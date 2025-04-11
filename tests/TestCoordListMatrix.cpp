#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>

#include "../include/CoordListMatrix.h"
#include "../include/Types.h"
#include <fstream>
#include <MatrixUtils.h>


// Helper function to create a small test matrix file
void createTestMatrixFileA(const std::string& filename) {
  std::ofstream fout(filename);
  if (!fout.is_open()) {
    throw std::runtime_error("Failed to open file A.");
  }

  fout << "2 3 3\n";

  fout << "1 2 3\n";
  fout << "2 2 2\n";
  fout << "2 3 1\n";
  fout.close();
}


// Helper function to create a small test matrix file
void createTestMatrixFileB(const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Failed to open file B.");
    }

    // M=3, N=7, NNZ=12
    fout << "3 7 12\n";

    // 1-based (row, col, val) lines
    fout << "1 2 4\n";
    fout << "2 2 1\n";
    fout << "1 3 3\n";
    fout << "2 3 1\n";
    fout << "3 3 1\n";
    fout << "2 4 1\n";
    fout << "3 4 2\n";
    fout << "2 5 2\n";
    fout << "3 5 2\n";
    fout << "2 6 2\n";
    fout << "3 6 1\n";
    fout << "3 7 1\n";

    fout.close();
}

TEST_CASE("CoordListMatrix constructor and getCoords", "[CoordListMatrix]") {
  const std::string filenameA = "Trec4.mtx";
  createTestMatrixFileA(filenameA);

  CoordListMatrix matrix(filenameA);

  std::pair<int, int> size = matrix.shape();
  int M = size.first;
  int N = size.second;

  // Check matrix dimensions
  REQUIRE(M == 2);
  REQUIRE(N == 3);

  // Expected coordinates
  // !!! 1-based converted to 0-based !!!
  std::vector<Coord> expectedCoords = {
    {0, 1},  // (1,2) -> (0,1)
    {1, 1},  // (2,2) -> (1,1)
    {1, 2}   // (2,3) -> (1,2)
  };

  // CHECK
  std::vector<Coord> actualCoords = matrix.getCoords();
  REQUIRE(actualCoords.size() == expectedCoords.size());

  // CHECK EVERY VALUE
  for (size_t i = 0; i < expectedCoords.size(); ++i) {
      REQUIRE(actualCoords[i] == expectedCoords[i]);
  }
}

TEST_CASE("CoordListMatrix constructor, vector", "[CoordListMatrix]") {
    int M = 1000;
    int N = 2000;

    SECTION("Negative dimensions throw invalid_argument") {
        std::vector<Coord> empty;
        REQUIRE_THROWS_AS(CoordListMatrix(empty, -1, N), std::invalid_argument);
        REQUIRE_THROWS_AS(CoordListMatrix(empty, M, -1), std::invalid_argument);
        REQUIRE_THROWS_AS(CoordListMatrix(empty, 0, 0), std::invalid_argument);
    }

    SECTION("Out-of-bounds coordinate throws out_of_range") {
        std::vector<Coord> badCoords = {
            {M, 10},    // row is out of bounds
            {100, N}, // col is out of bounds
            {-1, 10},   // negative row
            {10, -1}    // negative col
        };

        for (const auto& c : badCoords) {
            std::vector test = {c};
            REQUIRE_THROWS_AS(CoordListMatrix(test, M, N), std::out_of_range);
        }
    }

    SECTION("Valid coordinates construct successfully") {
        std::vector<Coord> goodCoords = {
            {0, 0},
            {M - 1, N - 1},
            {500, 1000}
        };

        REQUIRE_NOTHROW(CoordListMatrix(goodCoords, M, N));
    }
}

TEST_CASE("CoordListMatrix naive matmul dimensions match", "[CoordListMatrix]") {
    SECTION("Check mismatch error thrown") {
        double sparsity = 0.05;
        int M = 1000;
        int N = 2000;

        const auto coordsA = generateSparseMatrix(sparsity, M, 5);
        CoordListMatrix A(coordsA, M, 5);

        const auto coordsB = generateSparseMatrix(sparsity, 10, N);
        CoordListMatrix B(coordsB, 10, N);

        REQUIRE_THROWS_AS(A.naiveMatmul(B), std::invalid_argument);
    }

    SECTION("Check result dimensions") {
        double sparsity = 0.05;
        int M = 1000;
        int N = 2000;
        int K = 1000;

        const auto coordsA = generateSparseMatrix(sparsity, M, K);
        CoordListMatrix A(coordsA, M, K);

        const auto coordsB = generateSparseMatrix(sparsity, K, N);
        CoordListMatrix B(coordsB, K, N);

        auto result = A.naiveMatmul(B);
        REQUIRE(result.shape() == std::pair<int, int> (M, N));
    }
}


/**
 * @brief Test boolean matmul of Trec4 (2×3) and Trec5 (3×7).
 *        Result should be 2×7.
 */
TEST_CASE("CoordListMatrix matmul Trec4(2×3) vs Trec5(3×7)", "[CoordListMatrix]") {

    const std::string fileT4 = "Trec4.mtx";
    const std::string fileT5 = "Trec5.mtx";
    createTestMatrixFileA(fileT4);
    createTestMatrixFileB(fileT5);

    CoordListMatrix A(fileT4); // shape = 2×3
    CoordListMatrix B(fileT5); // shape = 3×7

    // Multiply A×B => shape 2×7
    CoordListMatrix C = A.naiveMatmul(B);
    auto [rowsC, colsC] = C.shape();

    REQUIRE(rowsC == 2);
    REQUIRE(colsC == 7);

    std::vector<Coord> expectedCoords = {
            // row=0
            {0,1},{0,2},{0,3},{0,4},{0,5},
            // row=1
            {1,1},{1,2},{1,3},{1,4},{1,5},{1,6}
    };

    std::vector<Coord> actualCoords = C.getCoords();

    // Sort them to match expectedCoords format
    auto coordSorter = [](const Coord &a, const Coord &b) {
        if (a.row != b.row) return a.row < b.row;
        return a.col < b.col;
    };
    std::sort(actualCoords.begin(), actualCoords.end(), coordSorter);
    std::sort(expectedCoords.begin(), expectedCoords.end(), coordSorter);

    REQUIRE(actualCoords.size() == expectedCoords.size());

    for (size_t i = 0; i < expectedCoords.size(); i++) {
        CHECK(actualCoords[i].row == expectedCoords[i].row);
        CHECK(actualCoords[i].col == expectedCoords[i].col);
    }
}