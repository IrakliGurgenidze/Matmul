#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>

#include "../include/CSRMatrix.h"
#include "../include/Types.h"
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
  std::vector<Coord> actualCoords;

  std::vector<int> rowPtr = csr.getRowPtr();
  std::vector<int> colIdx = csr.getColIdx();

  int rowPtrSize = static_cast<int>(rowPtr.size()) - 1;
  for (int row=0; row < rowPtrSize; ++row) {
    for (int i = rowPtr[row]; i < rowPtr[row + 1]; ++i) {
      int col = colIdx[i];
      actualCoords.push_back({row, col});
    }
  }


  REQUIRE(actualCoords.size() == expectedCoords.size());

  // CHECK EVERY VALUE
  for (size_t i = 0; i < expectedCoords.size(); ++i) {
    REQUIRE(actualCoords[i].row == expectedCoords[i].row);
    REQUIRE(actualCoords[i].col == expectedCoords[i].col);
  }
}


