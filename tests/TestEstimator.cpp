#include "../CoordListMatrix.h"
#include "../include/Estimator.h"
#include "../include/HashUtils.h"
#include "../include/MatrixUtils.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>

void benchmark_estimator(int M, int K, int N, double curr_sparsity,
                         const std::string &label, bool csv = false) {
  // std::random_device rd;
  int seedA = 1;
  int seedB = 2;

  const auto coordsA = generateSparseMatrix(curr_sparsity, M, K, seedA);
  const auto coordsB = generateSparseMatrix(curr_sparsity, K, N, seedB);

  CoordListMatrix A(coordsA, M, K);
  CoordListMatrix B(coordsB, K, N);

  const auto start = std::chrono::steady_clock::now();
  double estimation =
      estimateProductSize(A.getHashedCoords(), B.getHashedCoords());
  const auto end = std::chrono::steady_clock::now();

  const auto startGT = std::chrono::steady_clock::now();
  int trueNumberNonzero = groundTruthCalc(A.getCoords(), B.getCoords());
  const auto endGT = std::chrono::steady_clock::now();

  const auto elapsed = std::chrono::duration<double>(end - start).count();
  const auto elapsedGT = std::chrono::duration<double>(endGT - startGT).count();
  REQUIRE(elapsed > 0);

  if (!csv) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << std::left << std::setw(21) << label << std::right
              << "  Estimation=" << std::setw(15) << estimation
              << "  True NNZ=" << std::setw(15) << trueNumberNonzero
              << "  time=" << std::setw(10) << elapsed << " s" << std::endl;
  } else {
    const std::filesystem::path outDir{
        "/Users/griffinravo/CLionProjects/Matmul/output/sparsityest"};

    std::filesystem::create_directories(outDir);

    // Open output file in append mode
    std::ofstream fout(outDir / "SparsityEstvGT.csv", std::ios::app);
    if (!fout.is_open()) {
      throw std::runtime_error("Could not open CSV file");
    }

    fout << elapsed << "," << elapsedGT << "," << curr_sparsity << "\n";
    fout.close();
  }
}

TEST_CASE("Estimator returns k^2 when sketch cannot be filled", "[Estimator]") {
  // This test creates a small dataset with only 4 output (a, c) pairs:
  // - b = 1: (1, 1) joins with (1, 10), (1, 11) → (1, 10), (1, 11)
  // - b = 2: (2, 2), (3, 2) join with (2, 12) → (2, 12), (3, 12)
  // Total distinct (a, c) = 4

  std::vector<R1Tuple> R1 = {
      {1, 1, 0.1, 0.0}, // x=a=1, y=b=1, h1=0.1
      {2, 2, 0.2, 0.0}, // x=a=2, y=b=2
      {3, 2, 0.3, 0.0}, // x=a=3, y=b=2
  };
  std::vector<R2Tuple> R2 = {
      {1, 10, 0.0, 0.7}, // x=b=1, y=c=10, h2=0.7
      {1, 11, 0.0, 0.9}, // x=b=1, y=c=11
      {2, 12, 0.0, 0.6}, // x=b=2, y=c=12
  };

  double epsilon = 0.1;
  int k = static_cast<int>(9.0 / (epsilon * epsilon));
  double expectedUpperBound = k * k;

  double estimate = estimateProductSize(R1, R2, epsilon);

  REQUIRE(estimate == expectedUpperBound);
}

TEST_CASE("Estimator returns accurate value when sketch fills",
          "[Estimator][Integration]") {
  int seed = 42;
  double sparsity = 0.05;
  int n = 1000;

  std::vector<Coord> R1coords = generateSparseMatrix(sparsity, n, n, seed);
  std::vector<Coord> R2coords =
      generateSparseMatrix(sparsity, n, n, seed * 517);

  HashContext::instance().setSeeds(12345, 67890);

  std::vector<R1Tuple> R1;
  R1.reserve(R1coords.size());
  for (auto &[a, b] : R1coords) {
    R1.push_back({a, b, murmur_hash(a, HashContext::instance().seed1), 0.0});
  }

  std::vector<R2Tuple> R2;
  R2.reserve(R2coords.size());
  for (auto &[b, c] : R2coords) {
    R2.push_back({b, c, 0.0, murmur_hash(c, HashContext::instance().seed2)});
  }

  int groundTruth = groundTruthCalc(R1coords, R2coords);
  REQUIRE(groundTruth > 0);

  double epsilon = 0.1;
  int k = static_cast<int>(9.0 / (epsilon * epsilon));
  double estimate = estimateProductSize(R1, R2, epsilon);

  INFO("Ground truth size: " << groundTruth);
  INFO("Estimator result: " << estimate);
  INFO("k: " << k);

  if (groundTruth > k) {
    double lower = (1.0 - epsilon) * groundTruth;
    double upper = (1.0 + epsilon) * groundTruth;
    REQUIRE(estimate >= lower);
    REQUIRE(estimate <= upper);
  } else {
    REQUIRE(estimate <= static_cast<double>(k * k));
  }
}

TEST_CASE("Estimator run time", "[Estimator_SZ][SizeSweep]") {
  double sparsity = 0.00005;
  int start_N = 10000;
  int end_N = 250000;
  int step = 10000;

  for (int N = start_N; N <= end_N; N += step) {
    int M = N;
    int K = N; // or N, if you want square multiplications

    std::string label = "[estimator-sweep N=" + std::to_string(N) + "]";
    benchmark_estimator(M, K, N, sparsity, label, true);
  }
}

TEST_CASE("Estimator SP", "[Estimator_SP][Sparsity_Sweep]") {
  double start_sp = 0.00000000005;
  int N = 1000;

  for (double i = start_sp; i < 0.0000005; i *= 10) {
    int M = N;
    int K = N; // or N, if you want square multiplications
    std::string label = "[estimator-sweep N=" + std::to_string(N) +
                        " sp=" + std::to_string(i) + "]";
    benchmark_estimator(M, K, N, i, label, true);
  }
}
