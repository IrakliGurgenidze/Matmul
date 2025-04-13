#include "../include/CSRMatrix.h"
#include "../include/CoordListMatrix.h"
#include "../include/Estimator.h"

#include <chrono>
#include <iostream>
#include <vector>

int main() {
  initPairwiseHashes();

  // Load CoordList matrices (for estimation)
  CoordListMatrix matrix1("data/cage15.mtx");
  CoordListMatrix matrix2("data/cage15.mtx");

  auto R1Tuples = matrix1.getHashedCoords();
  auto R2Tuples = matrix2.getHashedCoords();

  // Estimate join-project size
  double estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);
  std::cout << "Estimated join-project size: " << estimate << std::endl;

  // Load CSR matrices
  CSRMatrix csrA("data/cage15.mtx");
  CSRMatrix csrB("data/cage15.mtx");

  // Time naive CSR matmul
  auto t1 = std::chrono::high_resolution_clock::now();
  CSRMatrix naiveResult = csrA.naiveMatmul(csrB);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> naiveTime = t2 - t1;
  std::cout << "Naive CSR matmul time: " << naiveTime.count() << " seconds"
            << std::endl;
  std::cout << "Naive result nnz: " << naiveResult.getCoords().size()
            << std::endl;

  // Time optimized CSR matmul (with estimate)
  auto t3 = std::chrono::high_resolution_clock::now();
  CSRMatrix optimizedResult = csrA.optimizedMatmul(csrB, estimate);
  auto t4 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> optimizedTime = t4 - t3;
  std::cout << "Optimized CSR matmul time: " << optimizedTime.count()
            << " seconds" << std::endl;
  std::cout << "Optimized result nnz: " << optimizedResult.getCoords().size()
            << std::endl;

  return 0;
}
