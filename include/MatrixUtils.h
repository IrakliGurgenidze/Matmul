#include <vector>
#include "Types.h"

#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

/*
 * Computes the ground-truth number of non-zero elements in the product of two matrices, given
 * their non-zero coordinates.
 *
 * @param r1: vector of non-zero Coords in first matrix
 * @param r2: vector of non-zero Coords in second matrix
 * @return: the number of non-zero elements in the product vector
 */
int groundTruthCalc(const std::vector<Coord>& r1, const std::vector<Coord>& r2);

/*
 * Generate a random sparse matrix based on provided parameters.
 *
 * @param sparseDegree: a ratio between 0 and 1 equaling non-zero / total elements in the result
 * @param numRows: number of rows in the result
 * @param numCols: number of cols in the result
 * @return: a vector of Coords representing the non-zero positions
 */
std::vector<Coord> generateSparseMatrix(double sparseDegree, int numRows, int numCols, int randomSeed=42);

#endif //MATRIXUTILS_H
