#include "../include/HashUtils.h"
#include "../include/external/MurmurHash3.h"
#include "../include/Types.h"
#include "../include/HashContext.h"
#include <random>
#include <cmath>



void initPairwiseHashes() {
    HashContext::instance().randomize();
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

double HashCoord::hAC() const  {return hashAC(h1,h2);}