#pragma once

#include <vector>
#include <string>
#include "Types.h"

#ifndef COORDLISTMATRIX_H
#define COORDLISTMATRIX_H


/**
 * @class CoordListMatrix
 * @brief Stores a sparse matrix as a list of non-zero (row, col) coordinates.
 *
 * Provides a simple coordinate list (COO) representation for use in sparse matrix operations.
 */
class CoordListMatrix{
  public:

    /**
     * @brief Loads non-zero entries from a Matrix Market (.mtx) file.
     *
     * Ignores zero entries and comment lines. Uses 0-based indexing.
     * Throws std::runtime_error on file or parse errors.
     */
    explicit CoordListMatrix(const std::string &filename);

    /**
     * @brief Returns the list of non-zero (row, col) coordinates.
     * @return Reference to the internal vector of Coord entries.
     */
    std::vector<Coord>& getCoords();

    /**
     * @brief Performs sparse matrix multiplication with this matrix on the left.
     *
     * Multiplies the current matrix (as left operand) with the given matrix `right`,
     * returning the coordinate list of the result.
     *
     * @param right The right-hand matrix in the multiplication (this × right)
     * @return CoordListMatrix representing the product
     *
     * Throws std::invalid_argument on matrix dimension mismatch.
     */
    CoordListMatrix matmul(CoordListMatrix &right);

    /**
     * @brief Returns the shape (numRows, numCols) of the matrix.
     * @return std::pair<int, int> representing (rows, cols)
     */
    [[nodiscard]] std::pair<int, int> shape() const;

  private:
    std::vector<Coord> coords;
    int M, N; // num rows, num cols
};

#endif //COORDLISTMATRIX_H
