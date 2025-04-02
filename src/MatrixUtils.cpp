#include <unordered_map>
#include <iostream>
#include "MatrixUtils.h"

void groundTruthCalc(std::vector<Coord> r1, std::vector<Coord> r2) {
    std::unordered_map<int, std::vector<int>> R1map;  // b -> all a's
    std::unordered_map<int, std::vector<int>> R2map;  // b -> all c's

    for (const auto &coord : r1) {
        R1map[coord.c].push_back(coord.r);  // R1: (a, b)
    }
    for (const auto &coord : r2) {
        R2map[coord.r].push_back(coord.c);  // R2: (b, c)
    }

    std::unordered_set<std::pair<int, int>, pair_hash> joinPairs;

    for (const auto &[b, aList] : R1map) {
        if (R2map.count(b)) {
            for (int a : aList) {
                for (int c : R2map[b]) {
                    joinPairs.emplace(a, c);
                }
            }
        }
    }

    std::cout << "Ground truth: " << joinPairs.size() << " join pairs\n";
}