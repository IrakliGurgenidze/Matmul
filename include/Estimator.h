#pragma once

#include "Types.h"
#include <vector>

#ifndef ESTIMATOR_H
#define ESTIMATOR_H

/**
 * @brief Replaces S with the smallest k values (by hash(ACpair)) in {S U F}.
 * Also clears F.
 *
 * @param S Vector containing smallest k ACpair (by hash value) discovered so
 * far
 * @param F Vector containing unprocessed ACpairs to be merged into S
 */
static void combine(std::vector<ACpair> &S, std::vector<ACpair> &F);

/**
 * @brief Performs a sweep-based merge step to update the top-k minimum hash
 * (a,c) pairs.
 *
 * This function scans the Cartesian product A × C for joinable tuples (where b
 * matches), computes their hashed (a,c) value using hAC = h1(a) - h2(c), and
 * updates:
 *   - S: the current set of the k smallest hash (a,c) pairs
 *   - F: the buffer of new candidate pairs that may enter S
 *
 * It avoids duplicate (a,c) pairs using the addedPairs set.
 * The threshold p is updated to the k-th smallest hash seen so far.
 *
 * @param A Vector of tuples from relation R1 (a, b, h1(a)).
 * @param C Vector of tuples from relation R2 (b, c, h2(c)).
 * @param p Hash threshold used to determine top-k inclusion (updated in-place).
 * @param k Number of minimum hash pairs to retain in S.
 * @param S Vector containing the smallest k ACpair entries (by hashAC)
 * discovered so far.
 * @param F Buffer of unprocessed ACpairs to be merged into S.
 * @param addedPairs Set used to avoid re-inserting duplicate (a, c) pairs.
 */
static void
pointerSweep(const std::vector<R1Tuple> &A, const std::vector<R2Tuple> &C,
             double &p, int k, std::vector<ACpair> &S, std::vector<ACpair> &F,
             std::unordered_set<std::pair<int, int>, pair_hash> &addedPairs);

/**
 * @brief Returns estimated number of non-zero values in product of sparse
 * matrix multiplication.
 *
 * @param R1in Vector of HashCoord, containing coordinates of non-zero values in
 * R1 and their prehashed values
 * @param R2in Vector of HashCoord, containing coordinates of non-zero values in
 * R2 and their prehashed values
 * @param epsilon Double in [0, 1) that sets error bound of the estimation
 *
 * @return Returns estimated number of non-zero values in product
 */
double estimateProductSize(const std::vector<HashCoord> &R1in,
                           const std::vector<HashCoord> &R2in,
                           double epsilon = 0.1);

#endif // ESTIMATOR_H
