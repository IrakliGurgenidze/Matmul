#pragma once

#include <utility>
#include <unordered_set>
#include <iostream>

#ifndef TYPES_H
#define TYPES_H

/**
 * @brief Represents an occupied coordinate in a sparse 2D matrix.
 */
struct Coord {
    int row, col;

    bool operator==(const Coord& other) const {
        return row == other.row && col == other.col;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Coord& c) {
    return os << "(" << c.row << ", " << c.col << ")";
}

/**
 * A tuple from relation R1 with a = row, b = join key, and pre-computed hash of 'a'.
 */
struct R1Tuple { int a, b; double h1a; };

/**
 * A tuple from relation R2 with b = join key, c = column, and precomputed hash of 'c'.
 */
struct R2Tuple { int b, c; double h2c; };

/**
 * Represents a pair (a, c) along with the hash value used in estimation.
 */
struct ACpair { int a, c; double hashAC; };

/**
 * Hash function for std::pair<int, int> to use in unordered_set/map.
 */
struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const;
};

#endif //TYPES_H
