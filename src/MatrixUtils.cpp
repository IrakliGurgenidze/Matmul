#include <unordered_map>
#include "MatrixUtils.h"

/*
 * Computes the ground-truth number of non-zero elements in the product of two matrices, given
 * their non-zero coordinates.
 *
 * @param r1: vector of non-zero Coords in first matrix
 * @param r2: vector of non-zero Coords in second matrix
 * @return: the number of non-zero elements in the product vector
 */
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