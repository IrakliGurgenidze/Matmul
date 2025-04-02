#include "../include/Estimator.h"
#include "../include/HashUtils.h"
#include <algorithm>
#include <iostream>
#include <unordered_map>

static double p_val = 1.0;
static int K_VAL = 0;  // Will be set in estimateProductSize

// Combines current sketch S with buffer F, keeping k smallest hashAC values
static void combine(std::vector<ACpair> &S, std::vector<ACpair> &F) {
    S.insert(S.end(), F.begin(), F.end());
    F.clear();

    if ((int)S.size() < K_VAL) return;

    std::nth_element(S.begin(), S.begin() + (K_VAL - 1), S.end(),
                     [](const ACpair &lhs, const ACpair &rhs) {
                         return lhs.hashAC < rhs.hashAC;
                     });

    double thresh = (S.begin() + (K_VAL - 1))->hashAC;

    auto mid = std::partition(S.begin(), S.end(), [&](const ACpair &x) {
        return x.hashAC <= thresh;
    });
    S.erase(mid, S.end());

    if ((int)S.size() > K_VAL) {
        S.resize(K_VAL);
    }

    p_val = thresh;
}

// Iterates over pairs from Ai x Ci and adds qualifying pairs to sketch
static void pointerSweep(const std::vector<R1Tuple> &A,
                         const std::vector<R2Tuple> &C,
                         double &p,
                         int k,
                         std::vector<ACpair> &S,
                         std::vector<ACpair> &F,
                         std::unordered_set<std::pair<int,int>, pair_hash> &addedPairs) {

    int s_bar = 0;
    for (int t = 0; t < (int)C.size(); t++) {
        while (hashAC(A[s_bar].h1a, C[t].h2c) > hashAC(A[(s_bar - 1 + A.size()) % A.size()].h1a, C[t].h2c)) {
            s_bar = (s_bar - 1 + A.size()) % A.size();
        }
        int s = s_bar;
        while (hashAC(A[s].h1a, C[t].h2c) < p) {
            double h = hashAC(A[s].h1a, C[t].h2c);
            std::pair<int,int> key(A[s].h1a, C[t].h2c);
            if (addedPairs.find(key) == addedPairs.end()) {
                addedPairs.insert(key);
                F.push_back({A[s].a, C[t].c, h});
            }
            if ((int)F.size() == k) {
                combine(S, F);
                p = p_val;
            }
            s = (s + 1) % A.size();
            if (s == s_bar) break;
        }
    }
}

// Main estimation function
// Assumes input tuples are already hashed
// Returns either k / p or an upper bound k^2 if sketch was not filled
double estimateProductSize(const std::vector<R1Tuple>& R1Tuples,
                           const std::vector<R2Tuple>& R2Tuples,
                           double epsilon) {
    initPairwiseHashes();
    K_VAL = static_cast<int>(9.0 / (epsilon * epsilon));
    p_val = 1.0;

    // Sort R1 and R2 by join key (b), then by hash value
    std::vector<R1Tuple> R1 = R1Tuples;
    std::vector<R2Tuple> R2 = R2Tuples;

    std::sort(R1.begin(), R1.end(), [](const R1Tuple &lhs, const R1Tuple &rhs) {
        return lhs.b != rhs.b ? lhs.b < rhs.b : lhs.h1a < rhs.h1a;
    });

    std::sort(R2.begin(), R2.end(), [](const R2Tuple &lhs, const R2Tuple &rhs) {
        return lhs.b != rhs.b ? lhs.b < rhs.b : lhs.h2c < rhs.h2c;
    });

    // Group R1 by b
    std::vector<std::pair<int, std::vector<R1Tuple>>> Ai;
    for (size_t i = 0; i < R1.size(); ) {
        int currB = R1[i].b;
        size_t j = i;
        while (j < R1.size() && R1[j].b == currB) j++;
        Ai.emplace_back(currB, std::vector<R1Tuple>(R1.begin() + i, R1.begin() + j));
        i = j;
    }

    // Group R2 by b
    std::vector<std::pair<int, std::vector<R2Tuple>>> Ci;
    for (size_t i = 0; i < R2.size(); ) {
        int currB = R2[i].b;
        size_t j = i;
        while (j < R2.size() && R2[j].b == currB) j++;
        Ci.emplace_back(currB, std::vector<R2Tuple>(R2.begin() + i, R2.begin() + j));
        i = j;
    }

    std::vector<ACpair> S, F;
    S.reserve(K_VAL);
    F.reserve(K_VAL);

    std::unordered_set<std::pair<int,int>, pair_hash> addedPairs;

    size_t i = 0, j = 0;
    while (i < Ai.size() && j < Ci.size()) {
        if (Ai[i].first == Ci[j].first) {
            pointerSweep(Ai[i].second, Ci[j].second, p_val, K_VAL, S, F, addedPairs);
            i++; j++;
        } else if (Ai[i].first < Ci[j].first) {
            i++;
        } else {
            j++;
        }
    }

    combine(S, F);

    if (static_cast<int>(S.size()) == K_VAL) {
        return static_cast<double>(K_VAL) / p_val;
    }
    return static_cast<double>(K_VAL) * K_VAL;
}
