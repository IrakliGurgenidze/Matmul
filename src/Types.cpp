#include "../include/Types.h"
#include <functional>

// Implementation of the hash function for std::pair<int, int>
std::size_t pair_hash::operator()(const std::pair<int, int> &p) const {
  return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
}
