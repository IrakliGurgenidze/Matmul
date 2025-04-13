#pragma once

#include "Types.h"
#include <vector>

#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

/**
 * @brief Computes the ground-truth number of non-zero elements in the product
 * of two matrices, given their non-zero coordinates.
 *
 * @param r1: Vector of non-zero Coords in first matrix
 * @param r2: Vector of non-zero Coords in second matrix
 *
 * @return: The number of non-zero elements in the product vector
 */
int groundTruthCalc(const std::vector<Coord> &r1, const std::vector<Coord> &r2);

/**
 * @brief Generates a random sparse matrix based on provided parameters.
 *
 * @param sparseDegree: A ratio between 0 and 1 equaling non-zero / total
 * elements in the result
 * @param numRows: Number of rows in the result
 * @param numCols: Number of cols in the result
 * @param randomSeed Seed for randomly generating the matrix
 *
 * @return: A vector of Coords representing the non-zero positions, sorted by
 * increasing row then col.
 */
std::vector<Coord> generateSparseMatrix(double sparseDegree, int numRows,
                                        int numCols, int randomSeed = 42);

#endif // MATRIXUTILS_H
