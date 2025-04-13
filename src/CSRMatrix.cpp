#include "../include/CSRMatrix.h"
#include <algorithm>
#include <Estimator.h>
#include <fstream>
#include <sstream>

// Comparator for sorting coords by row then col
static bool compareRowCol(const Coord &a, const Coord &b) {
  return a.row != b.row ? a.row < b.row : a.col < b.col;
}

static std::vector<int> visited;

inline void ensureVisitedSize( int numCols) {
  if (static_cast<int>(visited.size()) < numCols) {
    visited.assign(numCols, -1);
  }
}


CSRMatrix::CSRMatrix(const std::string &filename) {
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    throw std::runtime_error("Could not open " + filename);
  }

  std::string line;
  while (true) {
    if (!std::getline(fin, line)) {
      throw std::runtime_error("Could not read " + filename);
    }
    if (line.empty() || line[0] == '%') {
      continue;
    }

    std::stringstream ss(line);
    int nnz;
    ss >> this->M >> this->N >> nnz;
    if (ss.fail()) {
      throw std::runtime_error("Couldn't parse M, N, nnz from " + line);
    }

    // Read nnz lines of "row col [value]"
    std::vector<Coord> coords;
    coords.reserve(nnz);

    for (int i = 0; i < nnz; i++) {
      if (!std::getline(fin, line)) {
        throw std::runtime_error("File ended before reading all nonzeros.");
      }
      if (line.empty() || line[0] == '%') {
        i--;
        continue;
      }
      std::stringstream ssl(line);
      int r, c;
      double val;
      ssl >> r >> c >> val;
      if (ssl.fail()) {
        throw std::runtime_error("Couldn't parse a nonzero line: " + line);
      }

      // 1-based => 0-based
      r--;
      c--;

      if (val != 0.0) {
        coords.push_back({r, c});
      }
    }

    // CSR Conversion
    // Sort coords by row
    std::sort(coords.begin(), coords.end(), compareRowCol);

    // Set up row_ptr with 0's
    rowPtr.assign(this->M + 1, 0);

    // rowPtr setup
    for (const auto &[row, col] : coords) {
      rowPtr[row + 1]++;
    }
    for (int i = 0; i < M; ++i) {
      rowPtr[i + 1] += rowPtr[i];
    }

    // colIdx setup with offset from rowPtr
    colIdx.resize(coords.size());
    std::vector<int> rowOffset = rowPtr;
    for (const auto &[row, col] : coords) {
      int pos = rowOffset[row]++;
      colIdx[pos] = col;
    }
    return;
  }
}

CSRMatrix::CSRMatrix(const std::vector<Coord> &coords, int M, int N)
    : M(M), N(N) {

  if (M <= 0 || N <= 0) {
    throw std::invalid_argument("Matrix dimensions must be positive.");
  }

  for (const auto &[row, col] : coords) {
    if (row < 0 || row >= M || col < 0 || col >= N) {
      throw std::out_of_range("Coordinate is out of matrix bounds.");
    }
  }

  // Sort coords by (row, col) to ensure CSR canonical form
  std::vector<Coord> sortedCoords = coords;
  std::sort(sortedCoords.begin(), sortedCoords.end(),
            [](const Coord &a, const Coord &b) {
              return a.row != b.row ? a.row < b.row : a.col < b.col;
            });

  // Initialize rowPtr and colIdx
  rowPtr.assign(M + 1, 0);
  colIdx.resize(sortedCoords.size());

  // Count non-zeros per row
  for (const auto &[r, _] : sortedCoords) {
    rowPtr[r + 1]++;
  }

  // Exclusive prefix sum to get rowPtr
  for (int i = 0; i < M; ++i) {
    rowPtr[i + 1] += rowPtr[i];
  }

  // Fill colIdx
  std::vector<int> rowOffset = rowPtr; // Copy for placement tracking
  for (const auto &[r, c] : sortedCoords) {
    colIdx[rowOffset[r]++] = c;
  }
}

std::vector<Coord> CSRMatrix::getCoords() const {
  std::vector<Coord> coords;

  int rowPtrSize = static_cast<int>(rowPtr.size()) - 1;
  for (int row = 0; row < rowPtrSize; ++row) {
    for (int i = rowPtr[row]; i < rowPtr[row + 1]; ++i) {
      int col = colIdx[i];
      coords.push_back({row, col});
    }
  }

  return coords;
}

CSRMatrix CSRMatrix::naiveMatmul(const CSRMatrix &right) const {
  auto [rowsA, colsA] = this->shape();
  auto [rowsB, colsB] = right.shape();
  if (colsA != rowsB) {
    throw std::invalid_argument("matmul dimension mismatch: "
                                "Left cols (" +
                                std::to_string(colsA) + ") != Right rows (" +
                                std::to_string(rowsB) + ")");
  }

  // Row Pointers
  std::vector<int> resultRowPtr(rowsA + 1, 0);
  // Column Indices
  std::vector<int> resultColIdx;
  // Reserve space of last row val
  resultColIdx.reserve(resultRowPtr.back());

  // std::vector<std::unordered_set<int>> resSets(rowsA);
  ensureVisitedSize(colsB);

  for (int i = 0; i < rowsA; ++i) {
    const size_t before = resultColIdx.size();

    for (int aPos = rowPtr[i]; aPos < rowPtr[i + 1]; ++aPos) {
      // column index in A = row index in B
      int j = colIdx[aPos];

      for (int bPos = right.rowPtr[j]; bPos < right.rowPtr[j + 1]; ++bPos) {
        int k = right.colIdx[bPos];
        if (visited[k] != i) {
          visited[k] = i;
          resultColIdx.push_back(k);
        }
      }
    }
    std::sort(resultColIdx.begin() + before, resultColIdx.end());
    resultRowPtr[i + 1] = static_cast<int>(resultColIdx.size());
  }

  CSRMatrix result = *this;
  result.M = rowsA;
  result.N = colsB;
  result.rowPtr = std::move(resultRowPtr);
  result.colIdx = std::move(resultColIdx);

  return result;

}

CSRMatrix CSRMatrix::optimizedMatmul(const CSRMatrix &right, double estimate) {
  // Check for matrix dimension mismatch
  auto [rowsA, colsA] = this->shape();
  auto [rowsB, colsB] = right.shape();
  if (colsA != rowsB) {
    throw std::invalid_argument("matmul dimension mismatch: "
                                "Left cols (" +
                                std::to_string(colsA) + ") != Right rows (" +
                                std::to_string(rowsB) + ")");
  }

  // Row Pointers
  std::vector<int> resultRowPtr(rowsA + 1, 0);
  // Column Indices
  std::vector<int> resultColIdx;
  // Reserve space of last row val
  resultColIdx.reserve(static_cast<size_t>(estimate));

  // std::vector<std::unordered_set<int>> resSets(rowsA);
  ensureVisitedSize(colsB);


  for (int i = 0; i < rowsA; ++i) {
    const size_t before = resultColIdx.size();

    for (int aPos = rowPtr[i]; aPos < rowPtr[i + 1]; ++aPos) {
      // column index in A = row index in B
      int j = colIdx[aPos];

      for (int bPos = right.rowPtr[j]; bPos < right.rowPtr[j + 1]; ++bPos) {
        int k = right.colIdx[bPos];
        if (visited[k] != i) {
          visited[k] = i;
          resultColIdx.push_back(k);
        }
      }
    }
    std::sort(resultColIdx.begin() + before, resultColIdx.end());
    resultRowPtr[i + 1] = static_cast<int>(resultColIdx.size());
  }

  CSRMatrix result = *this;
  result.M = rowsA;
  result.N = colsB;
  result.rowPtr = std::move(resultRowPtr);
  result.colIdx = std::move(resultColIdx);

  return result;
}

std::vector<CSRMatrix> CSRMatrix::batchNaiveMatmul(
    const std::vector<CSRMatrix> &rights) const {
  auto [rowsA, colsA] = this->shape();

  // Validate all matrix dimensions first
  for (const auto &right : rights) {
    if (colsA != right.shape().first) {
      throw std::invalid_argument(
          "Dimension mismatch in naive batch multiplication.");
    }
  }

  std::vector<CSRMatrix> results;
  results.reserve(rights.size());

  // Now safe to run all multiplications
  for (const auto &right : rights) {
    results.push_back(this->naiveMatmul(right));
  }

  return results;
}



std::vector<CSRMatrix> CSRMatrix::batchOptimizedMatmul(
    const std::vector<CSRMatrix> &rights, double epsilon) const {

  for (const auto &right : rights) {
    if (this->N != right.shape().first) {
      throw std::invalid_argument("Dimension mismatch in batchOptimizedMatmul");
    }
  }

  std::vector<CSRMatrix> results;
  results.reserve(rights.size());

  std::vector<std::vector<int>> leftGroups(this->N);

  for (int row = 0; row < this->M; ++row) {
    for (int idx = rowPtr[row]; idx < rowPtr[row + 1]; ++idx) {
      int col = colIdx[idx];
      leftGroups[col].push_back(row);
    }
  }

  // Process each right matrix in the batch
  for (const auto &right : rights) {
    std::vector<std::vector<int>> rightGroups(right.M);
    for (int row = 0; row < right.M; ++row) {
      for (int idx = right.rowPtr[row]; idx < right.rowPtr[row+1]; ++idx) {
        int col = right.colIdx[idx];
        rightGroups[row].push_back(col);
      }
    }

    CoordListMatrix forEstimateA(this->getCoords(), this->M, this->N);
    auto [rightM, rightN] = right.shape();
    CoordListMatrix forEstimateB(right.getCoords(), rightM, rightN);

    // Call the estimator for the current left/right pair
    // The estimator uses the hashed coordinates from each matrix
    double estimatedJoinSize = estimateProductSize(
            forEstimateA.getHashedCoords(), forEstimateB.getHashedCoords(), epsilon);

    // Preallocate storage for the join result using the estimated join size.
    std::vector<Coord> resultCoords;
    resultCoords.reserve(static_cast<size_t>(estimatedJoinSize));

    // Instead of iterating by left row (as in naive), iterate over join keys
    // directly
    for (size_t joinKey = 0; joinKey < this->N; joinKey++) {
      if (!leftGroups[joinKey].empty() && !rightGroups[joinKey].empty()) {
        for (int leftRow : leftGroups[joinKey]) {
          for (int rightCol : rightGroups[joinKey]) {
            resultCoords.push_back({leftRow, rightCol});
          }
        }
      }
    }

    results.emplace_back(resultCoords, this->M, right.N);

  }
  return results;
}



std::pair<int, int> CSRMatrix::shape() const { return {this->M, this->N}; }
