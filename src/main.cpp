#include "../include/CoordListMatrix.h"
#include "../include/Estimator.h"

#include <iostream>
#include <vector>
#include <chrono>


int main() {
    // Initialize hash seeds prior to construction of matrices
    initPairwiseHashes();

    // Load sparse matrices from .mtx files
    CoordListMatrix matrix1("data/bwm200.mtx");
    CoordListMatrix matrix2("data/rdb200.mtx");

    std::vector<Coord> coords1 = matrix1.getCoords();
    std::vector<Coord> coords2 = matrix2.getCoords();


    auto R1Tuples = matrix1.getHashedCoords();
    auto R2Tuples = matrix2.getHashedCoords();

    // Estimate join-project size
    double estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);

    std::cout << "Estimated join-project size: " << estimate << std::endl;

    return 0;
}
