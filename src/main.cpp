#include "../include/CSRMatrix.h"
#include "../include/CoordListMatrix.h"
#include "../include/Estimator.h"
#include "../include/MatrixUtils.h"

#include <chrono>
#include <iostream>
#include <vector>
#include <filesystem>

int main() {
  // Initialize MurMur based pairwise hashes.
  initPairwiseHashes();


  // Demo with generated data
  std::cout << "--- DEMO USING GENERATED DATA ---" << std::endl;

  // Load CoordList matrices (for estimation)
  // Using two randomly generated datasets
  const auto coordsA = generateSparseMatrix(0.00005, 50000, 50000, 1);
  const auto coordsB = generateSparseMatrix(0.00005, 50000, 50000, 2);

  CoordListMatrix randMatrix1(coordsA, 50000, 50000);
  CoordListMatrix randMatrix2(coordsB, 50000, 50000);

  // Generate R1 and R2 sets
  auto R1Tuples = randMatrix1.getHashedCoords();
  auto R2Tuples = randMatrix2.getHashedCoords();

  // Estimator function to generate join-project estimate
  auto estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);
  std::cout << "Estimated join-project size: " << estimate << std::endl;

  // Timed CoordList matrix naive multiplication
  auto t1 = std::chrono::high_resolution_clock::now();
  auto resultCL = randMatrix1.naiveMatmul(randMatrix2);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = t2 - t1;

  std::cout << "CoordList naive matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCL.getCoords().size() << std::endl;

  // Timed CoordList matrix optimized multiplication
  t1 = std::chrono::high_resolution_clock::now();
  resultCL = randMatrix1.optimizedMatmul(randMatrix2, estimate);
  t2 = std::chrono::high_resolution_clock::now();
  elapsed = t2 - t1;

  std::cout << "CoordList optimized matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCL.getCoords().size() << std::endl;

  // Creating two CSR matrix representations from the two CoordLists
  CSRMatrix csrMatrix1 = CSRMatrix(randMatrix1.getCoords(), 50000, 50000);
  CSRMatrix csrMatrix2 = CSRMatrix(randMatrix2.getCoords(), 50000, 50000);

  // Timed CSR matrix naive multiplication
  t1 = std::chrono::high_resolution_clock::now();
  CSRMatrix resultCSR = csrMatrix1.naiveMatmul(csrMatrix2);
  t2 = std::chrono::high_resolution_clock::now();
  elapsed = t2 - t1;

  std::cout << "CSR naive matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCSR.getCoords().size() << std::endl;

  // Timed CSR matrix optimized multiplication
  t1 = std::chrono::high_resolution_clock::now();
  resultCSR = csrMatrix1.optimizedMatmul(csrMatrix2, estimate);
  t2 = std::chrono::high_resolution_clock::now();
  elapsed = t2 - t1;

  std::cout << "CSR optimized matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCSR.getCoords().size() << "\n" << std::endl;

  return 0;
}
