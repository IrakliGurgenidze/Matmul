#ifndef MATMUL_HASHCONTEXT_H
#define MATMUL_HASHCONTEXT_H

#pragma once
#include <cstdint>
#include <random>

inline constexpr uint64_t PRIME = 4294967311ULL;   // keep PRIME here

class HashContext {
public:
    uint64_t seed1{0}, seed2{0};

    static HashContext& instance() {
        static HashContext ctx;
        return ctx;
    }

    void randomize() {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> dis(1, PRIME - 1);
        seed1 = dis(gen);
        seed2 = dis(gen);
    }

    void setSeeds(uint64_t s1, uint64_t s2) { seed1 = s1; seed2 = s2; }
private:
    HashContext() = default;
};

#endif //MATMUL_HASHCONTEXT_H
