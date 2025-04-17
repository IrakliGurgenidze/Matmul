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

  std::cout << "--- DEMO USING REAL DATA ---" << std::endl;

  // Load CoordList matrices (for estimation)
  // Using two sample datasets
#include <filesystem>

    std::filesystem::path sourceDir = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path projectRoot = sourceDir.parent_path(); // go up from src/
    std::filesystem::path dataDir = projectRoot / "data";

    CoordListMatrix matrix1(dataDir / "bwm200.mtx");
    CoordListMatrix matrix2(dataDir / "rdb200.mtx");


  // Generate R1 and R2 sets
  auto R1Tuples = matrix1.getHashedCoords();
  auto R2Tuples = matrix2.getHashedCoords();

  // Estimator function to generate join-project estimate
  double estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);
  std::cout << "Estimated join-project size: " << estimate << std::endl;

  // Timed CoordList matrix naive multiplication
  auto t1 = std::chrono::high_resolution_clock::now();
  CoordListMatrix resultCL = matrix1.naiveMatmul(matrix2);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = t2 - t1;

  std::cout << "CoordList naive matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCL.getCoords().size() << std::endl;

  // Timed CoordList matrix optimized multiplication
  t1 = std::chrono::high_resolution_clock::now();
  resultCL = matrix1.optimizedMatmul(matrix2, estimate);
  t2 = std::chrono::high_resolution_clock::now();
  elapsed = t2 - t1;

  std::cout << "CoordList optimized matmul time: " << elapsed.count()
            << " seconds, nnz: " << resultCL.getCoords().size() << std::endl;

  // Creating two CSR matrix representations from the two CoordLists
  auto [M1, N1] = matrix1.shape();
  auto [M2, N2] = matrix2.shape();
  CSRMatrix csrMatrix1 = CSRMatrix(matrix1.getCoords(), M1, N1 );
  CSRMatrix csrMatrix2 = CSRMatrix(matrix2.getCoords(), M2, N2 );

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

  // *************************************************************

  // Demo with generated data
  std::cout << "--- DEMO USING GENERATED DATA ---" << std::endl;

  // Load CoordList matrices (for estimation)
  // Using two randomly generated datasets
  const auto coordsA = generateSparseMatrix(0.00005, 50000, 50000, 1);
  const auto coordsB = generateSparseMatrix(0.00005, 50000, 50000, 2);

  CoordListMatrix randMatrix1(coordsA, 50000, 50000);
  CoordListMatrix randMatrix2(coordsB, 50000, 50000);

  // Generate R1 and R2 sets
  R1Tuples = randMatrix1.getHashedCoords();
  R2Tuples = randMatrix2.getHashedCoords();

  // Estimator function to generate join-project estimate
  estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);
  std::cout << "Estimated join-project size: " << estimate << std::endl;

  // Timed CoordList matrix naive multiplication
  t1 = std::chrono::high_resolution_clock::now();
  resultCL = randMatrix1.naiveMatmul(randMatrix2);
  t2 = std::chrono::high_resolution_clock::now();
  elapsed = t2 - t1;

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
  csrMatrix1 = CSRMatrix(randMatrix1.getCoords(), 50000, 50000);
  csrMatrix2 = CSRMatrix(randMatrix2.getCoords(), 50000, 50000);

  // Timed CSR matrix naive multiplication
  t1 = std::chrono::high_resolution_clock::now();
  resultCSR = csrMatrix1.naiveMatmul(csrMatrix2);
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
