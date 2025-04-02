#include "../include/HashUtils.h"
#include "../include/Types.h"
#include "../include/external/MurmurHash3.h"
#include <random>
#include <cmath>

// Prime for generating random seeds
const uint64_t PRIME = 4294967311ULL;

// Global hash seeds
uint64_t murmurSeed1 = 0;
uint64_t murmurSeed2 = 0;

// Initializes murmur hash seeds with random values in [1, PRIME - 1]
void initPairwiseHashes() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, PRIME - 1);
    murmurSeed1 = dis(gen);
    murmurSeed2 = dis(gen);
}

// Applies MurmurHash3_x86_32 to integer x with given seed and maps to [0, 1]
double murmur_hash(int x, uint64_t seed) {
    uint32_t hash_val;
    MurmurHash3_x86_32(&x, sizeof(x), seed, &hash_val);
    return static_cast<double>(hash_val) / static_cast<double>(UINT32_MAX);
}

// Combines h1(a) and h2(c) into a single hash value in [0, 1]
double hashAC(double h1a, double h2c) {
    return std::fmod(h1a * 0.618 + h2c * 0.382, 1.0);
}