#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include "../include/Estimator.h"
#include "../include/MatrixUtils.h"
#include "../include/HashUtils.h"

TEST_CASE("Estimator returns k^2 when sketch cannot be filled", "[Estimator]") {
    // This test creates a small dataset with only 4 output (a, c) pairs:
    // - b = 1: (1, 1) joins with (1, 10), (1, 11) → (1, 10), (1, 11)
    // - b = 2: (2, 2), (3, 2) join with (2, 12) → (2, 12), (3, 12)
    // Total distinct (a, c) = 4

    std::vector<R1Tuple> R1 = {
        {1, 1, 0.1},  // a=1, b=1
        {2, 2, 0.2},  // a=2, b=2
        {3, 2, 0.3},  // a=3, b=2
    };

    std::vector<R2Tuple> R2 = {
        {1, 10, 0.7},  // b=1, c=10
        {1, 11, 0.9},  // b=1, c=11
        {2, 12, 0.6},  // b=2, c=12
    };

    double epsilon = 0.1;
    int k = static_cast<int>(9.0 / (epsilon * epsilon));
    double expectedUpperBound = k * k;

    double estimate = estimateProductSize(R1, R2, epsilon);

    REQUIRE(estimate == expectedUpperBound);
}


TEST_CASE("Estimator returns accurate value when sketch fills", "[Estimator][Integration]") {
    int seed = 42;
    double sparsity = 0.05;
    int n = 1000;

    std::vector<Coord> R1coords = generateSparseMatrix(sparsity, n, n, seed);
    std::vector<Coord> R2coords = generateSparseMatrix(sparsity, n, n, seed * 517);

    initPairwiseHashes();

    std::vector<R1Tuple> R1;
    for (const auto& [a, b] : R1coords) {
        R1.push_back({a, b, murmur_hash(a, murmurSeed1)});
    }

    std::vector<R2Tuple> R2;
    for (const auto& [b, c] : R2coords) {
        R2.push_back({b, c, murmur_hash(c, murmurSeed2)});
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