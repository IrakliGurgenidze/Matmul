#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include "../include/CoordListMatrix.h"
#include "../include/Types.h"
#include <Estimator.h>
#include <iomanip>
#include <filesystem>
#include <fstream>

TEST_CASE("File-based optimized matmul estimate vs ground truth", "[CoordListMatrix][Estimator][Integration]") {
    const std::string mtxFileA = "/Users/griffinravo/CLionProjects/Matmul/data/bwm200.mtx";
    const std::string mtxFileB = "/Users/griffinravo/CLionProjects/Matmul/data/rdb200.mtx";

    CoordListMatrix A(mtxFileA);
    CoordListMatrix B(mtxFileB);

    REQUIRE(A.shape().second == B.shape().first); // Ensure compatible dimensions

    double epsilon = 0.05;
    auto R1Tuples = A.getHashedCoords();
    auto R2Tuples = B.getHashedCoords();
    double estimatedNNZ = estimateProductSize(R1Tuples, R2Tuples, epsilon);

    CoordListMatrix C = A.optimizedMatmul(B, estimatedNNZ);
    int trueNNZ = static_cast<int>(C.getCoords().size());
    auto [rowsC, colsC] = C.shape();

    // Extract matrix names (strip .mtx extension)
    std::string nameA = std::filesystem::path(mtxFileA).stem().string();
    std::string nameB = std::filesystem::path(mtxFileB).stem().string();

    // CSV write
    const std::filesystem::path outDir{"/Users/griffinravo/CLionProjects/Matmul/output"};
    std::filesystem::create_directories(outDir);
    std::ofstream fout(outDir / "real_world_tests.csv", std::ios::app);
    REQUIRE(fout.is_open());

    fout << nameA << "," << nameB << "," << rowsC << "," << colsC << ","
         << trueNNZ << "," << estimatedNNZ << "," << epsilon << "\n";
    fout.close();

    INFO("rows = " << rowsC);
    INFO("cols = " << colsC);
    INFO("True NNZ = " << trueNNZ);
    INFO("Estimated NNZ = " << estimatedNNZ);

    REQUIRE(rowsC > 0);
    REQUIRE(colsC > 0);
    REQUIRE(trueNNZ >= 0);
    REQUIRE(estimatedNNZ >= 0.0);
}
