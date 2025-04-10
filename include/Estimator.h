#pragma once

#include <vector>
#include "Types.h"

#ifndef ESTIMATOR_H
#define ESTIMATOR_H

/**
 * @brief Replaces S with the smallest k values (by hash(ACpair)) in {S U F}. Also clears F.
 *
 * @param S Vector containing smallest k ACpair (by hash value) discovered so far
 * @param F Vector containing unprocessed ACpairs to be merged into S
 */
static void combine(std::vector<ACpair> &S, std::vector<ACpair> &F);

/**
 * @brief
 *
 * @param A
 * @param C
 * @param p
 * @param k
 * @param S
 * @param F
 * @param addedPairs
 */
static void pointerSweep(const std::vector<R1Tuple> &A,
                         const std::vector<R2Tuple> &C,
                         double &p,
                         int k,
                         std::vector<ACpair> &S,
                         std::vector<ACpair> &F,
                         std::unordered_set<std::pair<int,int>, pair_hash> &addedPairs);

/**
 * @brief Returns estimated number of non-zero values in product of sparse matrix multiplication.
 *
 * @param R1 Vector of R1Tuples, containing coordinates of non-zero values in R1 and their prehashed values
 * @param R2 Vector of R2Tuples, containing coordinates of non-zero values in R2 and their prehashed values
 * @param epsilon Double in [0, 1) that sets error bound of the estimation
 * @return Returns estimated number of non-zero values in product
 */
double estimateProductSize(std::vector<R1Tuple> &R1,
    std::vector<R2Tuple> &R2,
    double epsilon = 0.1);

#endif //ESTIMATOR_H
