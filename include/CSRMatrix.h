#ifndef CSRMATRIX_H
#define CSRMATRIX_H
#include <string>
#include <vector>


/**
 * @class CSRMatrix
 * @brief Stores a sparse matrix as 2 arrays,
 * where
 *
 * Provides a simple coordinate list (COO) representation for use in sparse matrix operations.
 */
class CSRMatrix {
  public:
    /**
     * @brief Loads non-zero entries from a Matrix Market (.mtx) file in CSR format.
     *
     * Ignores zero entries and comment lines. Uses 0-based indexing.
     * Throws std::runtime_error on file or parse errors.
     */
    explicit CSRMatrix(const std::string &filename);

    std::vector<int>& getRowPtr();
    std::vector<int>& getColIdx();

    /**
     * @brief Returns the shape (numRows, numCols) of the matrix.
     * @return std::pair<int, int> representing (rows, cols)
     */
    [[nodiscard]] std::pair<int, int> shape() const;



  private:
    std::vector<int> rowPtr;
    std::vector<int> colIdx;
    int M, N; // num rows, num cols
};



#endif //CSRMATRIX_H
