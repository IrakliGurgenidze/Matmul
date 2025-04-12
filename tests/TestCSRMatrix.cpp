#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>

#include "../include/CSRMatrix.h"
#include "../include/Types.h"
#include "../include/MatrixUtils.h"
#include <fstream>


// Helper function to create a small test matrix file
void CSRTestCreateTestMatrixFileA(const std::string& filename) {
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

// Helper function to create a small test matrix file
void CSRTestCreateTestMatrixFileB(const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Failed to open file B.");
    }

    // M=3, N=7, NNZ=12
    fout << "3 7 12\n";

    // 1-based (row, col, val) lines
    fout << "1 2 4\n";
    fout << "2 2 1\n";
    fout << "1 3 3\n";
    fout << "2 3 1\n";
    fout << "3 3 1\n";
    fout << "2 4 1\n";
    fout << "3 4 2\n";
    fout << "2 5 2\n";
    fout << "3 5 2\n";
    fout << "2 6 2\n";
    fout << "3 6 1\n";
    fout << "3 7 1\n";

    fout.close();
}


TEST_CASE("CSRMatrix constructor and getCSR", "[CSRMatrix]") {
    const std::string filenameA = "Trec4.mtx";
    CSRTestCreateTestMatrixFileA(filenameA);

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
    auto actualCoords = csr.getCoords();
    REQUIRE(actualCoords.size() == expectedCoords.size());

    // CHECK EVERY VALUE
    for (size_t i = 0; i < expectedCoords.size(); ++i) {
        REQUIRE(actualCoords[i] == expectedCoords[i]);
    }
}

TEST_CASE("CSRMatrix constructor, vector", "[CSRMatrix]") {
    int M = 1000;
    int N = 2000;

    SECTION("Negative dimensions throw invalid_argument") {
        std::vector<Coord> empty;
        REQUIRE_THROWS_AS(CSRMatrix(empty, -1, N), std::invalid_argument);
        REQUIRE_THROWS_AS(CSRMatrix(empty, M, -1), std::invalid_argument);
        REQUIRE_THROWS_AS(CSRMatrix(empty, 0, 0), std::invalid_argument);
    }

    SECTION("Out-of-bounds coordinate throws out_of_range") {
        std::vector<Coord> badCoords = {
            {M, 10},    // row is out of bounds
            {100, N}, // col is out of bounds
            {-1, 10},   // negative row
            {10, -1}    // negative col
        };

        for (const auto& c : badCoords) {
            std::vector test = {c};
            REQUIRE_THROWS_AS(CSRMatrix(test, M, N), std::out_of_range);
        }
    }

    SECTION("Valid coordinates construct successfully") {
        std::vector<Coord> goodCoords = {
            {0, 0},
            {M - 1, N - 1},
            {500, 1000}
        };

        REQUIRE_NOTHROW(CSRMatrix(goodCoords, M, N));
    }

    SECTION("Valid resulting coords") {
        double sparsity = 0.05;
        int M = 1000;
        int N = 2000;

        const auto actualCoords = generateSparseMatrix(sparsity, M, 5, 42);
        CSRMatrix mat(actualCoords, M, N);
        const auto matCoords = mat.getCoords();

        for (size_t i = 0; i < actualCoords.size(); ++i) {
            REQUIRE(actualCoords[i] == matCoords[i]);
        }
    }
}


TEST_CASE("CSRMatrix naiveMatmul dimensions match", "[CSRMatrix]") {
    SECTION("Check mismatch error thrown") {
        double sparsity = 0.05;
        int M = 1000;
        int N = 2000;

        const auto coordsA = generateSparseMatrix(sparsity, M, 5, 42);
        CSRMatrix A(coordsA, M, 5);

        const auto coordsB = generateSparseMatrix(sparsity, 10, N);
        CSRMatrix B(coordsB, 10, N);

        REQUIRE_THROWS_AS(A.naiveMatmul(B), std::invalid_argument);
    }

    SECTION("Check result dimensions") {
        double sparsity = 0.05;
        int M = 1000;
        int N = 2000;
        int K = 1000;

        const auto coordsA = generateSparseMatrix(sparsity, M, K);
        CSRMatrix A(coordsA, M, K);

        const auto coordsB = generateSparseMatrix(sparsity, K, N);
        CSRMatrix B(coordsB, K, N);

        auto result = A.naiveMatmul(B);
        REQUIRE(result.shape() == std::pair<int, int> (M, N));
    }
}



TEST_CASE("CSRMatrix naiveMatmul", "[CSRMatrix]") {
    const std::string fileT4 = "Trec4.mtx";
    const std::string fileT5 = "Trec5.mtx";
    CSRTestCreateTestMatrixFileA(fileT4);
    CSRTestCreateTestMatrixFileB(fileT5);

    CSRMatrix A(fileT4);
    CSRMatrix B(fileT5);

    // Multiply AxB => shape 2x7
    CSRMatrix C = A.naiveMatmul(B);
    auto [rowsC, colsC] = C.shape();

    REQUIRE(rowsC == 2);
    REQUIRE(colsC == 7);

    std::vector<Coord> expectedCoords = {
        // row=0
        {0,1},{0,2},{0,3},{0,4},{0,5},
        // row=1
        {1,1},{1,2},{1,3},{1,4},{1,5},{1,6}
    };

    std::vector<Coord> actualCoords = C.getCoords();

    // Sort them to match expectedCoords format
    auto coordSorter = [](const Coord &a, const Coord &b) {
        if (a.row != b.row) return a.row < b.row;
        return a.col < b.col;
    };
    std::sort(actualCoords.begin(), actualCoords.end(), coordSorter);
    std::sort(expectedCoords.begin(), expectedCoords.end(), coordSorter);

    REQUIRE(actualCoords.size() == expectedCoords.size());

    for (size_t i = 0; i < expectedCoords.size(); i++) {
        CHECK(actualCoords[i].row == expectedCoords[i].row);
        CHECK(actualCoords[i].col == expectedCoords[i].col);
    }
}
