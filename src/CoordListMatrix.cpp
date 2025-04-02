#include "../include/CoordListMatrix.h"
#include <fstream>
#include <sstream>

CoordListMatrix::CoordListMatrix(const std::string &filename, int &M, int &N) {
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
        ss >> M >> N >> nnz;
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
        this->coords = coords;
        return;
    }
}

std::vector<Coord>& CoordListMatrix::getCoords() {
    return coords;
}
