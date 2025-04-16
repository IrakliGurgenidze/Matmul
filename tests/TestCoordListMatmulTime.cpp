#define CATCH_CONFIG_MAIN // This tells Catch2 to provide a main() function
#include <CoordListMatrix.h>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <filesystem>
#include <fstream>

#include <Estimator.h>
#include <MatrixUtils.h>

void benchmarkCoordListMatmulNaive(int M, int K, int N, double sparsity,
                            const std::string &label) {
  int seedA = 1;
  int seedB = 2;

  const auto coordsA = generateSparseMatrix(sparsity, M, K, seedA);
  const auto coordsB = generateSparseMatrix(sparsity, K, N, seedB);

  CoordListMatrix A(coordsA, M, K);
  CoordListMatrix B(coordsB, K, N);

  const auto start = std::chrono::steady_clock::now();
  auto result = A.naiveMatmul(B);
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed = std::chrono::duration<double>(end - start).count();

  REQUIRE(elapsed > 0);

    const std::filesystem::path outDir{"/Users/griffinravo/CLionProjects/Matmul/output"};

    // Ensure output directory exists
    std::filesystem::create_directories(outDir);

    // Open output file in append mode
    std::ofstream fout(outDir / "CL_single_naive.csv", std::ios::app);
    if (!fout.is_open()) {
      throw std::runtime_error("Could not open CSV file");
    }

    fout << M << "," << K << "," << N << "," << elapsed << "\n";
    fout.close();

  REQUIRE(result.shape().first == M);
  REQUIRE(result.shape().second == N);
}

void benchmarkCoordListMatmulOptimized(int M, int K, int N, double sparsity,
                                const std::string &label) {
  // std::random_device rd;
  int seedA = 1;
  int seedB = 2;

  const auto coordsA = generateSparseMatrix(sparsity, M, K, seedA);
  const auto coordsB = generateSparseMatrix(sparsity, K, N, seedB);

  CoordListMatrix A(coordsA, M, K);
  CoordListMatrix B(coordsB, K, N);

  double estimation =
      estimateProductSize(A.getHashedCoords(), B.getHashedCoords());

  const auto start = std::chrono::steady_clock::now();
  auto result = A.optimizedMatmul(B, estimation);
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed = std::chrono::duration<double>(end - start).count();

  REQUIRE(elapsed > 0);

  const std::filesystem::path outDir{"/Users/griffinravo/CLionProjects/Matmul/output"};

  // Ensure output directory exists
  std::filesystem::create_directories(outDir);

  // Open output file in append mode
  std::ofstream fout(outDir / "CL_single_optimized.csv", std::ios::app);
  if (!fout.is_open()) {
    throw std::runtime_error("Could not open CSV file");
  }

  fout << M << "," << K << "," << N << "," << elapsed << "\n";
  fout.close();

  REQUIRE(result.shape().first == M);
  REQUIRE(result.shape().second == N);
}

void benchmarkCoordlistBatchedMatmulNaive(int N, double sparsity, int numMats, const std::string &label) {

  // Generate first matrix
  int seedA = 1;
  const auto coordsA = generateSparseMatrix(sparsity, N, N, seedA);
  CoordListMatrix A(coordsA, N, N);

  // Generate right-hand batch
  std::vector<CoordListMatrix> rights;
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

    const std::filesystem::path outDir{"/Users/griffinravo/CLionProjects/Matmul/output"};

    // Ensure output directory exists
    std::filesystem::create_directories(outDir);

    // Open output file in append mode
    std::ofstream fout(outDir / "CL_batched_naive.csv", std::ios::app);
    if (!fout.is_open()) {
      throw std::runtime_error("Could not open CSV file");
    }

    fout << N << "," << numMats << "," << elapsed << "\n";
    fout.close();

}

void benchmarkCoordlistBatchedMatmulOptimized(int N, double sparsity, int numMats, const std::string &label) {

  // Generate first matrix
  int seedA = 1;
  const auto coordsA = generateSparseMatrix(sparsity, N, N, seedA);
  CoordListMatrix A(coordsA, N, N);

  // Generate right-hand batch
  std::vector<CoordListMatrix> rights;
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

    const std::filesystem::path outDir{"/Users/griffinravo/CLionProjects/Matmul/output"};

    // Ensure output directory exists
    std::filesystem::create_directories(outDir);

    // Open output file in append mode
    std::ofstream fout(outDir / "CL_batched_optimized.csv", std::ios::app);
    if (!fout.is_open()) {
      throw std::runtime_error("Could not open CSV file");
    }

    fout << N << "," << numMats << "," << elapsed << "\n";
    fout.close();

}

// TEST_CASE("CoordList Performance Benchmark, (single op, naive)",
//           "[CoordListMatrix]") {
//   SECTION("Single matmul operation, small naive") {
//     benchmarkCoordListMatmulNaive(1000, 800, 1000, 0.005, "[naive matmul, small]");
//   }
//
//   SECTION("Single matmul operation, medium naive") {
//     benchmarkCoordListMatmulNaive(10000, 20000, 20000, 0.005,
//                            "[naive matmul, medium]");
//   }
//
//   SECTION("Single matmul operation, large naive") {
//     benchmarkCoordListMatmulNaive(300000, 300000, 3000000, 0.0005,
//                            "[naive matmul, large]");
//   }
// }

TEST_CASE("CoordList Scaling Sweep (single op, naive)", "[benchmark]") {
  double sparsity = 0.00005;
  int start_N = 20000;
  int end_N = 400000;
  int step = 20000;

  for (int N = start_N; N <= end_N; N += step) {
    int M = N;
    int K = N; // or N, if you want square multiplications

    std::string label = "[naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCoordListMatmulNaive(M, K, N, sparsity, label);
  }
}

TEST_CASE("CoordList Scaling Sweep (single optimized)", "[benchmark]") {
  double sparsity = 0.00005;
  int start_N = 20000;
  int end_N = 400000;
  int step = 20000;

  for (int N = start_N; N <= end_N; N += step) {
    int M = N;
    int K = N; // or N, if you want square multiplications

    std::string label = "[naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCoordListMatmulOptimized(M, K, N, sparsity, label);
  }
}

TEST_CASE("CoordList Scaling Sweep (batched naive)", "[benchmark][batch]") {
  double sparsity = 0.00005;
  int start_N = 20000;
  int end_N = 200000;
  int step = 20000;
  int numMats = 50; // Number of right-hand matrices in each batch

  for (int N = start_N; N <= end_N; N += step) {
    std::string label = "[batched-naive-sweep N=" + std::to_string(N) + "]";
    benchmarkCoordlistBatchedMatmulNaive(N, sparsity, numMats, label);
  }
}

TEST_CASE("CoordList Scaling Sweep (batched optimized)", "[benchmark][batch]") {
  double sparsity = 0.00005;
  int start_N = 20000;
  int end_N = 200000;
  int step = 20000;
  int numMats = 50; // Number of right-hand matrices in each batch

  for (int N = start_N; N <= end_N; N += step) {
    std::string label = "[batched-optimized-sweep N=" + std::to_string(N) + "]";
    benchmarkCoordlistBatchedMatmulOptimized(N, sparsity, numMats, label);
  }
}
