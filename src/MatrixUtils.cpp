#include <unordered_map>
#include <random>
#include <stdexcept>

#include "MatrixUtils.h"


int groundTruthCalc(const std::vector<Coord>& r1, const std::vector<Coord>& r2) {
    std::unordered_map<int, std::vector<int>> R1map;  // b -> all a's
    std::unordered_map<int, std::vector<int>> R2map;  // b -> all c's

    for (const auto &[r, c] : r1) R1map[c].push_back(r);  // A: (a, b)
    for (const auto &[r, c] : r2) R2map[r].push_back(c);  // B: (b, c)

    std::unordered_set<std::pair<int, int>, pair_hash> joinPairs;

    // perform join on b
    for (const auto &[b, aList] : R1map) {
        if (R2map.contains(b)) {
            for (int a : aList) {
                for (int c : R2map[b]) {
                    joinPairs.emplace(a, c);
                }
            }
        }
    }

    return static_cast<int>(joinPairs.size());
}

std::vector<Coord> generateSparseMatrix(double sparseDegree, int numRows, int numCols, int randomSeed) {
    if (sparseDegree <= 0.0 || sparseDegree > 1.0) {
        throw std::invalid_argument("sparseDegree must be in the range (0.0, 1.0]");
    }

    int totalElements = numRows * numCols;
    int targetNonzeros = static_cast<int>(totalElements * sparseDegree);

    std::mt19937 rng(randomSeed);
    std::uniform_int_distribution<int> rowDist(0, numRows - 1);
    std::uniform_int_distribution<int> colDist(0, numCols - 1);

    std::unordered_set<std::pair<int, int>, pair_hash> uniqueCoords;
    std::vector<Coord> result;

    while (static_cast<int>(uniqueCoords.size()) < targetNonzeros) {
        int r = rowDist(rng);
        int c = colDist(rng);
        std::pair<int, int> pos = {r, c};
        if (uniqueCoords.insert(pos).second) {
            result.push_back({r, c});
        }
    }

    return result;
}