#pragma once

#include <cstdint>
#include "HashContext.h"

#ifndef HASHUTILS_H
#define HASHUTILS_H

extern const uint64_t PRIME;

/**
 * @brief Initializes murmur hash seeds with random values in [1, PRIME - 1].
 */
void initPairwiseHashes();

/**
 * @brief Applies MurmurHash3_x86_32 to integer x with given seed and maps to [0, 1).
 * @param x Integer input to hash
 * @param seed A 64-bit seed for hash function randomization
 *
 * @return A double precision hash value uniformly distributed in [0, 1)
 */
double murmur_hash(int x, uint64_t seed);

/**
 * @brief Pairwise independent hash function for two pre-hashed values.
 *
 * @param h1a Hash value of a from R1
 * @param h2c Hash value of c from R2
 *
 * @return Hashed value of h1a and h2c, domain over [0,1)
 */
double hashAC(double h1a, double h2c);

#endif //HASHUTILS_H
