#include "../include/Estimator.h"
#include "../include/HashUtils.h"
#include <algorithm>
#include <iostream>
#include <unordered_map>

static double p_val = 1.0;
static int K_VAL = 0;  // Will be set in estimateProductSize

static void combine(std::vector<ACpair> &S, std::vector<ACpair> &F) {
    S.insert(S.end(), F.begin(), F.end());
    F.clear(); // Empty F

    if (static_cast<int>(S.size()) <= K_VAL) return;

    std::nth_element(S.begin(), S.begin() + (K_VAL - 1), S.end(),
                     [](const ACpair &lhs, const ACpair &rhs) {
                         return lhs.hAC() < rhs.hAC();
                     });

    double thresh = (S.begin() + (K_VAL - 1))->hAC();

    auto mid = std::partition(S.begin(), S.end(), [&](const ACpair &x) {
        return x.hAC() <= thresh;
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
                         std::unordered_set<uint64_t> &seen,
                         int maxC) {


    int Asize = A.size(); //cache a.size
    for (int t = 0; t < (int)C.size(); t++) { // loop through C
        double cHash = C[t].h2; //cache hash for col (yt)
        int s_bar = 0;
        // to find s_bar (starting row), find row in A w/ min hash value in this col (t)
        for (int s = 1; s < Asize; s++) {
            double currH = hashAC(A[s].h1, cHash); //hash (x_s_bar, y_t)
            double prevH = hashAC(A[s_bar].h1, cHash);
            if (currH < prevH) { //find s_bar s.t. hash is min
                s_bar  = s;
            }
        } // end PC line 12

        for (int i = 0; i < Asize; i++) {
            // Find all s where h(x,y) < p
            int s = (s_bar + i) % Asize; // line 19 check, gives s->(s_bar + offset) mod |A| where
            double h = hashAC(A[s].h1, cHash); // get hash
            if (h >= p) break; // only do while hash < p

            auto makeKey = [](int a, int c) {
                return (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(c);
            };
            // check for dups using seen vector, calculate proper idx
            uint64_t key = makeKey(A[s].x, C[t].y);
            if (seen.insert(key).second) {
                //push coord into F if not already
                F.push_back({A[s].x, C[t].y, A[s].h1, C[t].h2});
            }

            //When F reaches K capacity, combine, clear, update p locally
            if (F.size() >= K_VAL) {
                combine(S, F);
                p = p_val;
            }

        }
    }
}


double estimateProductSize(const std::vector<HashCoord> &R1in,
                           const std::vector<HashCoord> &R2in,
                           double epsilon) {


    K_VAL = static_cast<int>(9.0 / (epsilon * epsilon));
    std::vector<HashCoord> R1 = R1in;   // local copies
    std::vector<HashCoord> R2 = R2in;

    initPairwiseHashes();

    // Sort R1 and R2 by increasing join key (b), then by increasing hash value
    std::ranges::sort(R1, [](const R1Tuple &lhs, const R1Tuple &rhs) {
        return lhs.y != rhs.y ? lhs.y < rhs.y : lhs.h1 < rhs.h1;
    });

    std::ranges::sort(R2, [](const R2Tuple &lhs, const R2Tuple &rhs) {
        return lhs.x != rhs.x ? lhs.x < rhs.x : lhs.h2 < rhs.h2;
    });

    // Group R1 by b
    std::vector<std::pair<int, std::vector<R1Tuple>>> Ai;
    for (size_t i = 0; i < R1.size(); ) {
        int currB = R1[i].y;
        size_t j = i;
        while (j < R1.size() && R1[j].y == currB) j++;
        Ai.emplace_back(currB, std::vector(R1.begin() + i, R1.begin() + j));
        i = j;
    }

    // Group R2 by b
    std::vector<std::pair<int, std::vector<R2Tuple>>> Ci;
    for (size_t i = 0; i < R2.size(); ) {
        int currB = R2[i].x;
        size_t j = i;
        while (j < R2.size() && R2[j].x == currB) j++;
        Ci.emplace_back(currB, std::vector(R2.begin() + i, R2.begin() + j));
        i = j;
    }

    /*
    // Finding initial p value for O(n)
    // int maxProduct = 0;
    // for (size_t i = 0; i < Ai.size(); i++) {
    //     int product = static_cast<int>(Ai[i].second.size() * Ci[i].second.size());
    //     if (product > maxProduct) {
    //         maxProduct = product;
    //     }
    // }
    // p_val = std::min(1.0 / K_VAL, static_cast<double>(K_VAL) / maxProduct);
    */

    std::vector<HashCoord> S, F;
    S.reserve(K_VAL);
    F.reserve(K_VAL);

    // To check dupes, determine the max vals for a and c
    int maxA = 0, maxC = 0;
    for (auto& t : R1) maxA = std::max(maxA, t.x);
    for (auto& t : R2) maxC = std::max(maxC, t.y);
    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<std::size_t>(R1.size() + R2.size()));

    size_t i = 0, j = 0;
    while (i < Ai.size() && j < Ci.size()) {
        if (Ai[i].first == Ci[j].first) {
            pointerSweep(Ai[i].second, Ci[j].second, p_val, K_VAL, S, F, seen, maxC);
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
