#ifndef CSRMATRIX_H
#define CSRMATRIX_H

#include <CoordListMatrix.h>

#include "Types.h"
#include <string>
#include <vector>

/**
 * @class CSRMatrix
 * @brief Stores a sparse matrix in Compressed-Sparse-Row format.
 */
class CSRMatrix {
public:
  /**
   * @brief Loads non-zero entries from a Matrix Market (.mtx) file.
   *
   * Ignores zero entries and comment lines. Uses 0-based indexing.
   * @throws std::runtime_error on file or parse errors.
   */
  explicit CSRMatrix(const std::string &filename);

  /**
   * @brief Constructs a CSRMatrix from a vector of type Coord.
   *
   * @param coords Vector of Coord objects, representing occupied indices in the
   * matrix
   * @param M Number of rows in matrix
   * @param N Number of cols in matrix
   *
   * @throws std::invalid_argument on invalid M, N or mismatch between Coord
   * bounds and specified dimensions.
   */
  CSRMatrix(const std::vector<Coord> &coords, int M, int N);

  /**
   * @brief Returns a vector of non-zero (row, col) coordinates.
   * @return Vector of Coords listing non-zero indices in the matrix.
   */
  std::vector<Coord> getCoords() const;

  /**
   * @brief Performs naive (without product estimation) sparse matrix
   * multiplication with this matrix on the left.
   *
   * Multiplies the current matrix (as left operand) with the given matrix
   * `right`, returning the coordinate list of the result.
   *
   * @param right The right-hand matrix in the multiplication (this × right)
   * @return CSRMatrix representing the product
   *
   * @throws std::invalid_argument on matrix dimension mismatch.
   */
  CSRMatrix naiveMatmul(const CSRMatrix &right) const;

  /**
   * @brief Performs optimized (with estimation) sparse matrix multiplication
   * with this matrix on the left.
   *
   * Multiplies the current matrix (as left operand) with the given matrix
   * `right`, returning the coordinate list of the result.
   *
   * @param right The right-hand matrix in the multiplication (this × right)
   * @return CSRMatrix representing the product
   *
   * @throws std::invalid_argument on matrix dimension mismatch.
   */
  CSRMatrix optimizedMatmul(const CSRMatrix &right, double estimate);

  /**
   * @brief Performs naive batched sparse matrix multiplication
   * with this matrix on the left.
   *
   * Multiplies the current matrix (as left operand) with each matrix in
   * 'rights', returning a list of results.
   *
   * @param rights The right-hand matrices in the multiplication (this × right)
   * @return Vector of CSRMatrix representing the product
   *
   * @throws std::invalid_argument on matrix dimension mismatch.
   */
  [[nodiscard]] std::vector<CSRMatrix>
  batchNaiveMatmul(const std::vector<CSRMatrix> &rights) const;

  /**
   * @brief Performs optimized batched sparse matrix multiplication
   * with this matrix on the left.
   *
   * Multiplies the current matrix (as left operand) with each matrix in
   * 'rights', returning a list of results.
   *
   * @param rights The right-hand matrix in the multiplication (this × right)
   * @param epsilon The estimated product size of resulting matrix
   * @return Vector of CSRMatrix representing the product
   *
   * @throws std::invalid_argument on matrix dimension mismatch.
   */
  [[nodiscard]] std::vector<CSRMatrix>
  batchOptimizedMatmul(const std::vector<CSRMatrix> &rights,
                       double epsilon = 0.1) const;

  /**
   * @brief Returns the shape (numRows, numCols) of the matrix.
   * @return std::pair<int, int> representing (rows, cols)
   */
  [[nodiscard]] std::pair<int, int> shape() const;

private:
  // Only working with boolean matrices, so no values vector is needed.
  std::vector<int> rowPtr;
  std::vector<int> colIdx;

  int M, N; // num rows, num cols
};

#endif // CSRMATRIX_H
