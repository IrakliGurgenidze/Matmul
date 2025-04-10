#define CATCH_CONFIG_MAIN // This tells Catch2 to provide a main() function
#include <CoordListMatrix.h>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <chrono>

#include <MatrixUtils.h>

void benchmark_matmul_naive(
    int M, int K, int N,
    double sparsity,
    const std::string& label
) {
    const auto coordsA = generateSparseMatrix(sparsity, M, K);
    CoordListMatrix A(coordsA, M, K);

    const auto coordsB = generateSparseMatrix(sparsity, K, N);
    CoordListMatrix B(coordsB, K, N);

    const auto start = std::chrono::high_resolution_clock::now();
    auto result = A.naiveMatmul(B);
    const auto end = std::chrono::high_resolution_clock::now();

    const auto elapsed = std::chrono::duration<double>(end - start).count();

    const auto [rowsA, colsA] = A.shape();
    const auto [rowsB, colsB] = B.shape();

    std::cout << "RESULTS " + label + "\n";
    std::cout << "A shape: (" << rowsA << ", " << colsA << ")\n";
    std::cout << "B shape: (" << rowsB << ", " << colsB << ")\n";
    std::cout << "Elapsed time: " << elapsed << " seconds\n\n";

    REQUIRE(result.shape().first == M);  // sanity check on result dims
    REQUIRE(result.shape().second == N);
}

TEST_CASE("CoordList Performance Benchmark, (single op, naive)", "[CoordListMatrix]") {
    SECTION("Single matmul operation, small naive") {
        benchmark_matmul_naive(1000, 800, 1000, 0.005, "[naive matmul, small]");
    }

    SECTION("Single matmul operation, medium naive") {
        benchmark_matmul_naive(10000, 20000, 20000, 0.005, "[naive matmul, medium]");
    }

    SECTION("Single matmul operation, large naive") {
        benchmark_matmul_naive(100000, 50000, 100000, 0.005, "[naive matmul, medium]");
    }

}