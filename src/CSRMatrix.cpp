#include "../include/CSRMatrix.h"
#include <fstream>
#include <sstream>
#include <algorithm>

#include "Types.h"

bool compareRowCol(const Coord& a, const Coord& b) {
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


std::pair<int, int> CSRMatrix::shape() const {
    return {this->M, this->N};
}

std::vector<int>& CSRMatrix::getRowPtr() {
    return rowPtr;
}
std::vector<int>& CSRMatrix::getColIdx() {
    return colIdx;
}