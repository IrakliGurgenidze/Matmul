#include "../include/HashUtils.h"
#include "../include/external/MurmurHash3.h"
#include <random>
#include <cmath>

// Prime for generating random seeds
constexpr uint64_t PRIME = 4294967311ULL;

// Global hash seeds
uint64_t murmurSeed1 = 0;
uint64_t murmurSeed2 = 0;

void initPairwiseHashes() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, PRIME - 1);
    murmurSeed1 = dis(gen);
    murmurSeed2 = dis(gen);
}

double murmur_hash(int x, uint64_t seed) {
    uint32_t hash_val;
    MurmurHash3_x86_32(&x, sizeof(x), seed, &hash_val);
    return static_cast<double>(hash_val) / static_cast<double>(UINT32_MAX);
}

double hashAC(double h1a, double h2c) {
    double diff = h1a - h2c;
    if (diff < 0) diff += 1.0;
    return diff;
}