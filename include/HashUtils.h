#pragma once

#include <cstdint>

#ifndef HASHUTILS_H
#define HASHUTILS_H

extern const uint64_t PRIME;

// Initializes random seeds for hash functions used in estimation
void initPairwiseHashes();

// Hashes an integer x using MurmurHash3 with a given seed and maps to [0, 1]
double murmur_hash(int x, uint64_t seed);

// Combines two hash values into a final hash in [0, 1] for (a, c) pair
double hashAC(double h1a, double h2c);

// Global hash seeds shared across modules (initialized via initPairwiseHashes)
extern uint64_t murmurSeed1;
extern uint64_t murmurSeed2;

#endif //HASHUTILS_H
