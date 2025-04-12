#include "../include/CoordListMatrix.h"
#include "../include/Estimator.h"
#include <fstream>
#include <sstream>

CoordListMatrix::CoordListMatrix(const std::string &filename) : M(0), N(0){
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("Could not open " + filename);
    }

    // Skip comment lines
    std::string line;
    while (true) {
        if (!std::getline(fin, line)) {
            throw std::runtime_error("Invalid file: no dimension line found.");
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
        hashedCoords_.reserve(nnz);
        auto& ctx = HashContext::instance();

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
                hashedCoords_.push_back({r,c, murmur_hash(r, ctx.seed1), murmur_hash(c, ctx.seed2)});
            }
        }
        this->coords = coords;
        return;
    }
}

CoordListMatrix::CoordListMatrix(const std::vector<Coord>& coords, int M, int N)
    : M(M), N(N), coords(coords)
{
    if (M <= 0 || N <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive.");
    }

    for (const auto&[row, col] : coords) {
        if (row < 0 || row >= M || col < 0 || col >= N) {
            throw std::out_of_range("Coordinate is out of matrix bounds.");
        }
    }

    auto& ctx = HashContext::instance();
    hashedCoords_.reserve(coords.size());
    for (auto& p : coords) {
        hashedCoords_.push_back({
                                        p.row, p.col,
                                        murmur_hash(p.row, ctx.seed1),
                                        murmur_hash(p.col, ctx.seed2)
                                });
    }
}

std::vector<Coord>& CoordListMatrix::getCoords() {
    return coords;
}

const std::vector<HashCoord>& CoordListMatrix::getHashedCoords() const { return hashedCoords_; }


CoordListMatrix CoordListMatrix::naiveMatmul(const CoordListMatrix &right) {
    if (this->N != right.M) {
        throw std::invalid_argument("Matrix dimensions do not align.");
    }

    // Check for matrix dimension mismatch
    auto [rowsA, colsA] = this->shape();
    auto [rowsB, colsB] = right.shape();
    if (colsA != rowsB) {
        throw std::invalid_argument("matmul dimension mismatch: "
                                    "Left cols (" + std::to_string(colsA) + ") != Right rows (" + std::to_string(rowsB) + ")");
    }

    // Build row-based adjacency
    std::vector<std::vector<int>> rowA(rowsA);
    std::vector<std::vector<int>> rowB(rowsB);

    for(auto &c: coords){
        rowA[c.row].push_back(c.col);
    }

    for(auto &c: right.coords){
        rowB[c.row].push_back(c.col);
    }

    std::vector<std::unordered_set<int>> resSets(rowsA);

    for(int i = 0; i < rowA.size(); i++) {
        for (int j: rowA[i]) {
            if (j >= 0 && j < rowsB) {
                for (int k: rowB[j]) {
                    resSets[i].insert(k);
                }
            }
        }
    }

    std::vector<Coord> resultCoords;
    resultCoords.reserve(colsA);
    for (int i = 0; i < rowsA; i++) {
        for (int k : resSets[i]) {
            resultCoords.push_back({i, k});
        }
    }

    CoordListMatrix outMatrix = *this;
    outMatrix.coords = std::move(resultCoords);
    outMatrix.M = rowsA;
    outMatrix.N = colsB;

    return outMatrix;
}

CoordListMatrix CoordListMatrix::optimizedMatmul(const CoordListMatrix &right, double estimate) {
    if (this->N != right.M) {
        throw std::invalid_argument("Matrix dimensions do not align.");
    }

    // Check for matrix dimension mismatch
    auto [rowsA, colsA] = this->shape();
    auto [rowsB, colsB] = right.shape();
    if (colsA != rowsB) {
        throw std::invalid_argument("matmul dimension mismatch: "
                                    "Left cols (" + std::to_string(colsA) + ") != Right rows (" + std::to_string(rowsB) + ")");
    }

    // Preallocate the result coordinate vector using the externally provided estimate
    std::vector<Coord> resultCoords;
    resultCoords.reserve(static_cast<size_t>(estimate));

    // Build row-based adjacency
    std::vector<std::vector<int>> rowA(rowsA);
    for (const auto &c : this->coords) {
        rowA[c.row].push_back(c.col);
    }

    std::vector<std::vector<int>> rowB(rowsB);
    for (const auto &c : right.coords) {
        rowB[c.row].push_back(c.col);
    }

    std::vector<std::unordered_set<int>> resSets(rowsA);
    for (int i = 0; i < rowsA; i++) {
        for (int j : rowA[i]) {
            if (j >= 0 && j < rowsB) {
                for (int k : rowB[j]) {
                    resSets[i].insert(k);
                }
            }
        }
    }


    for (int i = 0; i < rowsA; i++) {
        for (int k : resSets[i]) {
            resultCoords.push_back({i, k});
        }
    }

    // Build and return the output CoordListMatrix
    CoordListMatrix outMatrix = *this;
    outMatrix.coords = std::move(resultCoords);
    outMatrix.M = rowsA;
    outMatrix.N = colsB;

    return outMatrix;
}

std::pair<int, int> CoordListMatrix::shape() const {
    return {this->M, this->N};
}
