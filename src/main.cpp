#include "../include/CoordListMatrix.h"
#include "../include/Estimator.h"
#include "../include/HashUtils.h"

#include <iostream>
#include <vector>

int main() {
    int M = 200, N = 200;

    // Load sparse matrices from .mtx files
    // Load sparse matrices from .mtx files
    CoordListMatrix matrix1("data/bwm200.mtx");
    CoordListMatrix matrix2("data/rdb200.mtx");

    std::vector<Coord> coords1 = matrix1.getCoords();
    std::vector<Coord> coords2 = matrix2.getCoords();

    // Initialize hash seeds
    initPairwiseHashes();

    // Convert Coord -> R1Tuple
    std::vector<R1Tuple> R1Tuples;
    for (const auto& coord : coords1) {
        R1Tuples.push_back({coord.row, coord.col, murmur_hash(coord.row, murmurSeed1)});
    }

    // Convert Coord -> R2Tuple
    std::vector<R2Tuple> R2Tuples;
    for (const auto& coord : coords2) {
        R2Tuples.push_back({coord.row, coord.col, murmur_hash(coord.col, murmurSeed2)});
    }

    // Estimate join-project size
    double estimate = estimateProductSize(R1Tuples, R2Tuples, 0.1);

    std::cout << "Estimated join-project size: " << estimate << std::endl;

    return 0;
}
