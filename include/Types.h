#pragma once

#include "HashUtils.h"
#include <iostream>
#include <unordered_set>
#include <utility>

#ifndef TYPES_H
#define TYPES_H

/**
 * @brief Represents an occupied coordinate in a sparse 2D matrix.
 */
struct Coord {
  int row, col;

  bool operator==(const Coord &other) const {
    return row == other.row && col == other.col;
  }
};

/**
 * Unified struct containing row/col values and their hashes
 */
struct HashCoord {
  int row;
  int col;
  double h1;
  double h2;

  double hAC() const;
};

using R1Tuple = HashCoord;
using R2Tuple = HashCoord;
using ACpair = HashCoord;

inline std::ostream &operator<<(std::ostream &os, const Coord &c) {
  return os << "(" << c.row << ", " << c.col << ")";
}

/**
 * Hash function for std::pair<int, int> to use in unordered_set/map.
 */
struct pair_hash {
  std::size_t operator()(const std::pair<int, int> &p) const;
};

#endif // TYPES_H
