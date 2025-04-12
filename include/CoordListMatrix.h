#ifndef COORDLISTMATRIX_H
#define COORDLISTMATRIX_H

#include <vector>
#include <string>
#include "Types.h"

/**
 * @class CoordListMatrix
 * @brief Stores a sparse matrix as a vector of non-zero (row, col) coordinates.
 */
class CoordListMatrix{
  public:

    /**
     * @brief Loads non-zero entries from a Matrix Market (.mtx) file.
     *
     * Ignores zero entries and comment lines. Uses 0-based indexing.
     * @throws std::runtime_error on file or parse errors.
     */
    explicit CoordListMatrix(const std::string &filename);

    /**
     * @brief Constructs a CoordList matrix from a vector of type Coord.
     *
     * @param coords Vector of Coord objects, representing occupied indices in the matrix
     * @param M Number of rows in matrix
     * @param N Number of cols in matrix
     *
     * @throws std::invalid_argument on invalid M, N or mismatch between Coord bounds and specified dimensions.
     */
    CoordListMatrix(const std::vector<Coord> &coords, int M, int N);

    /**
     * @brief Returns a vector of non-zero (row, col) coordinates.
     * @return Reference to the internal vector of Coord entries.
     */
    std::vector<Coord>& getCoords();

    /**
     * @brief Performs naive (without product estimation) sparse matrix multiplication with this matrix on the left.
     *
     * Multiplies the current matrix (as left operand) with the given matrix `right`,
     * returning the coordinate list of the result.
     *
     * @param right The right-hand matrix in the multiplication (this × right)
     * @return CoordListMatrix representing the product
     *
     * @throws std::invalid_argument on matrix dimension mismatch.
     */
    CoordListMatrix naiveMatmul(const CoordListMatrix &right);

    /**
     * @brief Performs optimized (with estimation) sparse matrix multiplication with this matrix on the left.
     *
     * Multiplies the current matrix (as left operand) with the given matrix `right`,
     * returning the coordinate list of the result.
     *
     * @param right The right-hand matrix in the multiplication (this × right)
     * @param estimation The estimated product size of resulting matrix
     * @return CoordListMatrix representing the product
     *
     * @throws std::invalid_argument on matrix dimension mismatch.
     */
    CoordListMatrix optimizedMatmul(const CoordListMatrix &right, double estimation);


    /**
     * @brief Performs optimized (with estimation) sparse matrix multiplication with this matrix on the left.
     *
     * Multiplies the current matrix (as left operand) with the given matrix `right`,
     * returning the coordinate list of the result.
     *
     * @param right The right-hand matrix in the multiplication (this × right)
     * @param estimation The estimated product size of resulting matrix
     * @return CoordListMatrix representing the product
     *
     * @throws std::invalid_argument on matrix dimension mismatch.
     */
      [[nodiscard]] std::vector<CoordListMatrix> batchOptimizedMatmul(const std::vector<CoordListMatrix>& rights,
        double epsilon = 0.1) const;


        /**
       * @brief Performs optimized (with estimation) sparse matrix multiplication with this matrix on the left.
       *
       * Multiplies the current matrix (as left operand) with the given matrix `right`,
       * returning the coordinate list of the result.
       *
       * @param right The right-hand matrix in the multiplication (this × right)
       * @param estimation The estimated product size of resulting matrix
       * @return CoordListMatrix representing the product
       *
       * @throws std::invalid_argument on matrix dimension mismatch.
       */
      [[nodiscard]] std::vector<CoordListMatrix> batchedNaiveMatmul(const std::vector<CoordListMatrix>& rights) const;


    /**
    * @brief Returns a vector of non-zero (row, col) coordinates with their respective hash values.
    * @return Reference to the internal vector of HashCoord entries.
    */
    const std::vector<HashCoord>& getHashedCoords() const;


    /**
     * @brief Returns the shape (numRows, numCols) of the matrix.
     * @return std::pair<int, int> representing (rows, cols)
     */
    [[nodiscard]] std::pair<int, int> shape() const;


  private:
    std::vector<Coord> coords;
    std::vector<HashCoord> hashedCoords_;
    int M, N; // num rows, num cols
};

#endif //COORDLISTMATRIX_H
