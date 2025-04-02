#pragma once

#include <utility>
#include <unordered_set>

#ifndef TYPES_H
#define TYPES_H

// Represents a nonzero coordinate in a sparse matrix (row, column)
struct Coord { int r, c; };

// A tuple from relation R1 with a = row, b = join key, and precomputed hash of 'a'
struct R1Tuple { int a, b; double h1a; };

// A tuple from relation R2 with b = join key, c = column, and precomputed hash of 'c'
struct R2Tuple { int b, c; double h2c; };

// Represents a pair (a, c) along with the hash value used in estimation
struct ACpair { int a, c; double hashAC; };

// Hash function for std::pair<int, int> to use in unordered_set/map
struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const;
};

#endif //TYPES_H
