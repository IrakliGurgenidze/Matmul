#include "../include/CSRMatrix.h"
#include <fstream>
#include <sstream>
#include <algorithm>

// Comparator for sorting coords by row then col
static bool compareRowCol(const Coord& a, const Coord& b) {
    return a.row != b.row ? a.row < b.row : a.col < b.col;
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
        for (const auto&[row, col] : coords) {
            rowPtr[row + 1]++;
        }
        for (int i = 0; i < M; ++i) {
            rowPtr[i + 1] += rowPtr[i];
        }

        // colIdx setup with offset from rowPtr
        colIdx.resize(coords.size());
        std::vector<int> rowOffset = rowPtr;
        for (const auto&[row, col] : coords) {
            int pos = rowOffset[row]++;
            colIdx[pos] = col;
        }
        return;
    }
}

CSRMatrix::CSRMatrix(const std::vector<Coord>& coords, int M, int N)
    : M(M), N(N) {

    if (M <= 0 || N <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive.");
    }

    for (const auto&[row, col] : coords) {
        if (row < 0 || row >= M || col < 0 || col >= N) {
            throw std::out_of_range("Coordinate is out of matrix bounds.");
        }
    }

    // Sort coords by (row, col) to ensure CSR canonical form
    std::vector<Coord> sortedCoords = coords;
    std::sort(sortedCoords.begin(), sortedCoords.end(),
        [](const Coord& a, const Coord& b) {
            return a.row != b.row ? a.row < b.row : a.col < b.col;
        });

    // Initialize rowPtr and colIdx
    rowPtr.assign(M + 1, 0);
    colIdx.resize(sortedCoords.size());

    // Count non-zeros per row
    for (const auto& [r, _] : sortedCoords) {
        rowPtr[r + 1]++;
    }

    // Exclusive prefix sum to get rowPtr
    for (int i = 0; i < M; ++i) {
        rowPtr[i + 1] += rowPtr[i];
    }

    // Fill colIdx
    std::vector<int> rowOffset = rowPtr;  // Copy for placement tracking
    for (const auto& [r, c] : sortedCoords) {
        colIdx[rowOffset[r]++] = c;
    }
}

std::vector<Coord> CSRMatrix::getCoords() const {
    std::vector<Coord> coords;

    int rowPtrSize = static_cast<int>(rowPtr.size()) - 1;
    for (int row=0; row < rowPtrSize; ++row) {
        for (int i = rowPtr[row]; i < rowPtr[row + 1]; ++i) {
            int col = colIdx[i];
            coords.push_back({row, col});
        }
    }

    return coords;
}


CSRMatrix CSRMatrix::naiveMatmul(const CSRMatrix &right) {
    // Check for matrix dimension mismatch
    auto [rowsA, colsA] = this->shape();
    auto [rowsB, colsB] = right.shape();
    if (colsA != rowsB) {
        throw std::invalid_argument("matmul dimension mismatch: "
                                    "Left cols (" + std::to_string(colsA) + ") != Right rows (" + std::to_string(rowsB) + ")");
    }

    std::vector<int> resultRowPtr(rowsA + 1, 0);
    std::vector<std::unordered_set<int>> resSets(rowsA);

    for (int i = 0; i < rowsA; ++i) {
        for (int aPos = rowPtr[i]; aPos < rowPtr[i + 1]; ++aPos) {
            // column index in A = row index in B
            int j = colIdx[aPos];

            for (int bPos = right.rowPtr[j]; bPos < right.rowPtr[j + 1]; ++bPos) {
                int k = right.colIdx[bPos];
                resSets[i].insert(k);
            }
        }
        resultRowPtr[i + 1] = static_cast<int>(resSets[i].size());
    }

    // prefix sum to finalize resultRowPtr
    for (int i = 0; i < rowsA; ++i) {
        resultRowPtr[i + 1] += resultRowPtr[i];
    }

    std::vector<int> resultColIdx;
    // Reserve space of last row val
    resultColIdx.reserve(resultRowPtr.back());

    // flatten resSets into colIdx
    for (int i = 0; i < rowsA; ++i) {
        for (int k: resSets[i]) {
            resultColIdx.push_back(k);
        }
    }

    CSRMatrix result = *this;
    result.M = rowsA;
    result.N = colsB;
    result.rowPtr = resultRowPtr;
    result.colIdx = resultColIdx;

    return result;
    // return CSRMatrix(std::vector<Coord>{}, this->M, right.N);
}

std::pair<int, int> CSRMatrix::shape() const {
    return {this->M, this->N};
}
