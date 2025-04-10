#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>

#include "../include/CSRMatrix.h"
#include "../include/Types.h"
#include "../include/MatrixUtils.h"
#include <fstream>


// Helper function to create a small test matrix file
void createTestMatrixFile(const std::string& filename) {
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


TEST_CASE("CSRMatrix constructor and getCSR", "[CSRMatrix]") {
  const std::string filenameA = "Trec4.mtx";
  createTestMatrixFile(filenameA);

  CSRMatrix csr(filenameA);

  std::pair<int, int> size = csr.shape();
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
  auto actualCoords = csr.getCoords();
  REQUIRE(actualCoords.size() == expectedCoords.size());

  // CHECK EVERY VALUE
  for (size_t i = 0; i < expectedCoords.size(); ++i) {
    REQUIRE(actualCoords[i] == expectedCoords[i]);
  }
}

TEST_CASE("CSRMatrix constructor, vector", "[CSRMatrix]") {
  int M = 1000;
  int N = 2000;

  SECTION("Negative dimensions throw invalid_argument") {
    std::vector<Coord> empty;
    REQUIRE_THROWS_AS(CSRMatrix(empty, -1, N), std::invalid_argument);
    REQUIRE_THROWS_AS(CSRMatrix(empty, M, -1), std::invalid_argument);
    REQUIRE_THROWS_AS(CSRMatrix(empty, 0, 0), std::invalid_argument);
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
      REQUIRE_THROWS_AS(CSRMatrix(test, M, N), std::out_of_range);
    }
  }

  SECTION("Valid coordinates construct successfully") {
    std::vector<Coord> goodCoords = {
      {0, 0},
      {M - 1, N - 1},
      {500, 1000}
    };

    REQUIRE_NOTHROW(CSRMatrix(goodCoords, M, N));
  }

  SECTION("Valid resulting coords") {
    double sparsity = 0.05;
    int M = 1000;
    int N = 2000;

    const auto actualCoords = generateSparseMatrix(sparsity, M, 5);
    CSRMatrix mat(actualCoords, M, N);
    const auto matCoords = mat.getCoords();

    for (size_t i = 0; i < actualCoords.size(); ++i) {
      REQUIRE(actualCoords[i] == matCoords[i]);
    }
  }
}
