#define CATCH_CONFIG_MAIN // This tells Catch2 to provide a main() function
#include <CSRMatrix.h>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdint>
#include <iomanip>

#include <Estimator.h>
#include <MatrixUtils.h>

void benchmarkCSRMatmulNaive(int M, int K, int N, double sparsity,
                            const std::string &label, bool csv = false) {
  int seedA = 1;
  int seedB = 2;

  const auto coordsA = generateSparseMatrix(sparsity, M, K, seedA);
  const auto coordsB = generateSparseMatrix(sparsity, K, N, seedB);

  CSRMatrix A(coordsA, M, K);
  CSRMatrix B(coordsB, K, N);

  const auto start = std::chrono::steady_clock::now();
  auto result = A.naiveMatmul(B);
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed = std::chrono::duration<double>(end - start).count();

  REQUIRE(elapsed > 0);

  if (!csv) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(21) << label << std::right
              << "  M=" << std::setw(7) << M << "  K=" << std::setw(7) << K
              << "  N=" << std::setw(7) << N << "  time=" << std::setw(10)
              << elapsed << " s" << std::endl;
  } else {
    std::cout << M << "," << K << "," << N << "," << elapsed << std::endl;
  }

  REQUIRE(result.shape().first == M);
  REQUIRE(result.shape().second == N);
}

void benchmarkCSRMatmulOptimized(int M, int K, int N, double sparsity,
                                const std::string &label, bool csv = false) {
  // std::random_device rd;
  int seedA = 1;
  int seedB = 2;

  const auto coordsA = generateSparseMatrix(sparsity, M, K, seedA);
  const auto coordsB = generateSparseMatrix(sparsity, K, N, seedB);

  CSRMatrix A(coordsA, M, K);
  CSRMatrix B(coordsB, K, N);

  CoordListMatrix forEstimateA(coordsA, M, K);
  CoordListMatrix forEstimateB(coordsB, K, N);

  double estimation =
      estimateProductSize(forEstimateA.getHashedCoords(), forEstimateB.getHashedCoords());

  const auto start = std::chrono::steady_clock::now();
  auto result = A.optimizedMatmul(B, estimation);
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed = std::chrono::duration<double>(end - start).count();

  REQUIRE(elapsed > 0);

  if (!csv) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(21) << label << std::right
              << "  M=" << std::setw(7) << M << "  K=" << std::setw(7) << K
              << "  N=" << std::setw(7) << N << "  time=" << std::setw(10)
              << elapsed << " s" << std::endl;
  } else {
    std::cout << M << "," << K << "," << N << "," << elapsed << std::endl;
  }

  // Measure result against ground truth, and verify dimensions match.
  // int trueNumberNonzero = groundTruthCalc(A.getCoords(), B.getCoords());
  // int computedNumberNonzero = result.getCoords().size();
  // REQUIRE(computedNumberNonzero == trueNumberNonzero);

  REQUIRE(result.shape().first == M);
  REQUIRE(result.shape().second == N);
}

void benchmarkCSRBatchedMatmulNaive(int N, double sparsity, int numMats, const std::string &label,
  bool csv = false) {

  // Generate first matrix
  int seedA = 1;
  const auto coordsA = generateSparseMatrix(sparsity, N, N, seedA);
  CSRMatrix A(coordsA, N, N);

  // Generate right-hand batch
  std::vector<CSRMatrix> rights;
  for (int i = 0; i < numMats; ++i) {
    int seedB = 100 + i;
    auto coordsB = generateSparseMatrix(sparsity, N, N, seedB);
    rights.emplace_back(coordsB, N, N);
  }

  // Time the batched naive matmul
  const auto start = std::chrono::steady_clock::now();
  auto results = A.batchNaiveMatmul(rights);
  const auto end = std::chrono::steady_clock::now();
  const double elapsed = std::chrono::duration<double>(end - start).count();

  // Verify results
  for (const auto& result : results) {
    REQUIRE(result.shape() == std::pair<int, int>{N, N});
  }

  // Print results
  if (csv) {
    std::cout << N << "," << numMats << "," << elapsed << std::endl;
  } else {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left
              << std::setw(25) << label
              << std::right
              << "  N=" << std::setw(6) << N
              << "  B=" << std::setw(4) << numMats
              << "  time=" << std::setw(10) << elapsed << " s"
              << std::endl;
  }
}

void benchmarkCSRBatchedMatmulOptimized(int N, double sparsity, int numMats, const std::string &label,
  bool csv = false) {

  // Generate first matrix
  int seedA = 1;
  const auto coordsA = generateSparseMatrix(sparsity, N, N, seedA);
  CSRMatrix A(coordsA, N, N);

  // Generate right-hand batch
  std::vector<CSRMatrix> rights;
  for (int i = 0; i < numMats; ++i) {
    int seedB = 100 + i;
    auto coordsB = generateSparseMatrix(sparsity, N, N, seedB);
    rights.emplace_back(coordsB, N, N);
  }

  // Time the batched optimized matmul
  const auto start = std::chrono::steady_clock::now();
  auto results = A.batchOptimizedMatmul(rights);
  const auto end = std::chrono::steady_clock::now();
  const double elapsed = std::chrono::duration<double>(end - start).count();

  // Verify results
  for (const auto& result : results) {
    REQUIRE(result.shape() == std::pair<int, int>{N, N});
  }

  // Print results
  if (csv) {
    std::cout << N << "," << numMats << "," << elapsed << std::endl;
  } else {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left
              << std::setw(25) << label
              << std::right
              << "  N=" << std::setw(6) << N
              << "  B=" << std::setw(4) << numMats
              << "  time=" << std::setw(10) << elapsed << " s"
              << std::endl;
  }
}

/** CSR matmul TESTS **/

TEST_CASE("CSR Performance Benchmark, (single op, naive)",
          "[CSRMatrix]") {
  SECTION("Single matmul operation, small naive") {
    benchmarkCSRMatmulNaive(1000, 800, 1000, 0.005, "[naive matmul, small]");
  }

  SECTION("Single matmul operation, medium naive") {
    benchmarkCSRMatmulNaive(10000, 20000, 20000, 0.005,
                           "[naive matmul, medium]");
  }

  // SECTION("Single matmul operation, large naive") {
  //   benchmarkCSRMatmulNaive(300000, 300000, 3000000, 0.0005,
  //                          "[naive matmul, large]");
  // }
}

TEST_CASE("CSR Scaling Sweep (single op, naive)", "[benchmark]") {
  // double sparsity = 0.001;
  double sparsity = 0.0005;

  int start_N = 10000;
  int end_N = 100000;
  int step = 10000;

  // int start_N = 1000;
  // int end_N = 10000;
  // int step = 1000;

  for (int N = start_N; N <= end_N; N += step) {
    int M = N;
    int K = N; // or N, if you want square multiplications

    std::string label = "[naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCSRMatmulNaive(M, K, N, sparsity, label, true);
  }
}

TEST_CASE("CSR Scaling Sweep (single op, optimized)", "[benchmark]") {
  // double sparsity = 0.001;
  double sparsity = 0.0005;

  int start_N = 10000;
  int end_N = 100000;
  int step = 10000;

  // int start_` 1000;

  for (int N = start_N; N <= end_N; N += step) {
    int M = N;
    int K = N; // or N, if you want square multiplications

    std::string label = "[naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCSRMatmulOptimized(M, K, N, sparsity, label, true);
  }
}



TEST_CASE("CSR Scaling Sweep (batched naive)", "[benchmark][batch]") {
  double sparsity = 0.00005;
  int start_N = 10000;
  int end_N = 100000;
  int step = 10000;
  int numMats = 20; // Number of right-hand matrices in each batch


  // int start_N = 1000;
  // int end_N = 10000;
  // int step = 1000;

  for (int N = start_N; N <= end_N; N += step) {
    std::string label = "[batched-naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCSRBatchedMatmulNaive(N, sparsity, numMats, label, true);
  }
}

TEST_CASE("CSR Scaling Sweep (batched optimized)", "[benchmark][batch]") {
  double sparsity = 0.00005;
  int start_N = 10000;
  int end_N = 100000;
  int step = 10000;
  int numMats = 20; // Number of right-hand matrices in each batch

  // int start_N = 1000;
  // int end_N = 10000;
  // int step = 1000;

  for (int N = start_N; N <= end_N; N += step) {
    std::string label = "[batched-optimized-sweep N=" + std::to_string(N) + "]";
    benchmarkCSRBatchedMatmulOptimized(N, sparsity, numMats, label, true);
  }
}
