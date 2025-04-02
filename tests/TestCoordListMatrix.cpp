#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>

#include "../include/CoordListMatrix.h"
#include "../include/Types.h"
#include <fstream>


// Helper function to create a small test matrix file
void createTestMatrixFile(const std::string& filename) {
  std::ofstream fout(filename);
  if (!fout.is_open()) {
    throw std::runtime_error("Failed to open file.");
  }
  fout << "2 3 3\n";  // Dimensions: 2x3 matrix with 3 nonzero values
  fout << "1 2 3\n";
  fout << "2 2 2\n";
  fout << "2 3 1\n";
  fout.close();

  //2 3 3
  //1 2 3
  //2 2 2
  //2 3 1

}

TEST_CASE("CoordListMatrix constructor and getCoords", "[CoordListMatrix]") {
  const std::string filename = "Trec4.mtx";
  createTestMatrixFile(filename);

  int M, N;

  CoordListMatrix matrix(filename, M, N);

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
    REQUIRE(actualCoords[i].row == expectedCoords[i].row);
    REQUIRE(actualCoords[i].col == expectedCoords[i].col);
  }
}
