#define CATCH_CONFIG_MAIN // This tells Catch2 to provide a main() function
#include <CoordListMatrix.h>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <chrono>
#include <iomanip>

#include <MatrixUtils.h>

void benchmark_matmul_naive(
    int M, int K, int N,
    double sparsity,
    const std::string& label,
    bool csv=false
) {
    const auto coordsA = generateSparseMatrix(sparsity, M, K);
    CoordListMatrix A(coordsA, M, K);

    const auto coordsB = generateSparseMatrix(sparsity, K, N);
    CoordListMatrix B(coordsB, K, N);

    const auto start = std::chrono::high_resolution_clock::now();
    auto result = A.naiveMatmul(B);
    const auto end = std::chrono::high_resolution_clock::now();

    const auto elapsed = std::chrono::duration<double>(end - start).count();

    if (!csv) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::left
                  << std::setw(21) << label
                  << std::right
                  << "  M=" << std::setw(7) << M
                  << "  K=" << std::setw(7) << K
                  << "  N=" << std::setw(7) << N
                  << "  time=" << std::setw(10) << elapsed << " s"
                  << std::endl;
    } else {
        std::cout << M << "," << K << "," << N << "," << elapsed << std::endl;
    }

    REQUIRE(result.shape().first == M);
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
        benchmark_matmul_naive(100000, 50000, 100000, 0.005, "[naive matmul, large]");
    }

}

TEST_CASE("CoordList Scaling Sweep (single op, naive)", "[benchmark]") {
    double sparsity = 0.0005;
    int start_N = 5000;
    int end_N = 750000;
    int step = 5000;

    for (int N = start_N; N <= end_N; N += step) {
        int M = N;
        int K = N / 1; // or N, if you want square multiplications

        std::string label = "[naive-sweep N=" + std::to_string(N) + "]";
        benchmark_matmul_naive(M, K, N, sparsity, label, true);
    }
}
