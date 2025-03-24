#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <unordered_set>

#include "MurmurHash3.h"

double EPSILON = 0.1;
int K_VAL = (int) (9.0 / (EPSILON*EPSILON));
const uint64_t PRIME = 4294967311ULL;

double p_val = 1.0;

struct Coord {
    int r;
    int c;
};

struct R1Tuple {
  int a;
  int b;
  double h1a;
};

struct R2Tuple {
  int b;
  int c;
  double h2c;
};

struct ACpair {
    int a;
    int c;
    double hashAC;
};

struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
    }
};

std::vector<Coord> readMtxAsCoordList(const std::string &filename, int &M, int &N) {
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
        return coords;
    }
}

uint64_t murmurSeed1, murmurSeed2;

void initPairwiseHashes(){
    std::random_device rd;
    std::mt19937_64 gen((rd()));
    std::uniform_int_distribution<uint64_t> dis(1, PRIME-1);
    murmurSeed1 = dis(gen);
    murmurSeed2 = dis(gen);
}

double murmur_hash(int x, uint64_t seed) {
    uint32_t hash_val;
    MurmurHash3_x86_32(&x, sizeof(x), seed, &hash_val);
    return static_cast<double>(hash_val) / static_cast<double>(UINT32_MAX);
}


double hashAC(double h1a, double h2c){
    return std::fmod(h1a * 0.618 + h2c * 0.382, 1.0);
}


void combine(std::vector<ACpair> &S, std::vector<ACpair> &F) {

    // h(S) union h(F)
    for(auto &x: F) S.push_back(x);
    F.clear();

    // No kth smallest val if size < k
    if(S.size() < K_VAL) return;


    // "RANK": split S with nth element func, lhs < k (at proper idx) < rhs
    // psuedocode: gives us 'v'
    std::nth_element(S.begin(), S.begin() + (K_VAL-1),S.end(),
                     [](const ACpair &lhs, const ACpair &rhs) {
        return lhs.hashAC < rhs.hashAC;
    });

    double thresh = (S.begin() + (K_VAL-1)) -> hashAC;

    // x in S where h(x) <= v
    auto mid = std::partition(S.begin(), S.end(),
                         [&](auto &x){return x.hashAC <= thresh;}
    );
    S.erase(mid, S.end());

    // if there are more than k items <= newp, keep exactly k
    if ((int)S.size() > K_VAL) {
        S.resize(K_VAL);
    }

    p_val = thresh;
}

void pointerSweep(
    const std::vector<R1Tuple> &A,
    const std::vector<R2Tuple> &C,
    double &p,              // threshold
    int k,
    std::vector<ACpair> &S, // global set of k smallest so far
    std::vector<ACpair> &F  // buffer that fills up to k
){


  int s_bar = 0;
  for(int t = 0; t < C.size(); t++) {
     while(hashAC(A[s_bar].h1a, C[t].h2c) > hashAC(A[(s_bar-1) % A.size()].h1a, C[t].h2c)){
       s_bar = ((s_bar-1) + A.size()) % A.size();
     }
     int s = s_bar;
     while(hashAC(A[s].h1a, C[t].h2c) < p){
       double h1 = hashAC(A[s].h1a, C[t].h2c);
       F.push_back({A[s].a, C[t].c, h1});
       if(F.size() == K_VAL) {
         combine(S, F);

         p = p_val;
       }
       s = (s + 1) % A.size();
       if(s == s_bar) break;
     }
  }
}

void groundTruthCalc(std::vector<Coord> r1, std::vector<Coord> r2) {
    std::unordered_map<int, std::vector<int>> R1map;  // b -> all a's
    std::unordered_map<int, std::vector<int>> R2map;  // b -> all c's

    for (const auto &coord : r1) {
        R1map[coord.c].push_back(coord.r);  // R1: (a, b)
    }
    for (const auto &coord : r2) {
        R2map[coord.r].push_back(coord.c);  // R2: (b, c)
    }

    std::unordered_set<std::pair<int, int>, pair_hash> joinPairs;

    for (const auto &[b, aList] : R1map) {
        if (R2map.count(b)) {
            for (int a : aList) {
                for (int c : R2map[b]) {
                    joinPairs.emplace(a, c);
                }
            }
        }
    }

    std::cout << "Ground truth: " << joinPairs.size() << " join pairs\n";

}



int main() {
    int M = 200;
    int N = 200;

    std::vector<Coord> R1 = readMtxAsCoordList("/Users/griffinravo/CLionProjects/Matmul/data/bwm200.mtx", M, N);
    std::vector<Coord> R2 = readMtxAsCoordList("/Users/griffinravo/CLionProjects/Matmul/data/rdb200.mtx", M, N);

    groundTruthCalc(R1, R2);

    std::vector<R1Tuple> R1Tuples;
    std::vector<R2Tuple> R2Tuples;

    initPairwiseHashes();

    // R1 -> (a,b,h1(a))
    for (auto &coord : R1) {
        int a = coord.r;
        int b = coord.c;
        double hval = murmur_hash(a, murmurSeed1);
        R1Tuples.push_back({a, b, hval});
    }

    //R2 -> (b,c,h2(c))
    for (auto &coord : R2) {
        int b = coord.r;
        int c = coord.c;
        double hval = murmur_hash(c, murmurSeed2);
        R2Tuples.push_back({b, c, hval});
    }


    // Sorting primarily by b then by hash val
    sort(R1Tuples.begin(), R1Tuples.end(), [](const R1Tuple &lhs, const R1Tuple &rhs) {
      if(lhs.b != rhs.b) return lhs.b < rhs.b;
      return lhs.h1a < rhs.h1a;
    });

    sort(R2Tuples.begin(), R2Tuples.end(), [](const R2Tuple &lhs, const R2Tuple &rhs) {
        if(lhs.b != rhs.b) return lhs.b < rhs.b;
        return lhs.h2c < rhs.h2c;
    });

    //Grouping by b for all groups of b=i
    std::vector<std::pair<int,std::vector<R1Tuple>>> Ai;
    int start = 0;

    while(start < R1Tuples.size()){
        int currB = R1Tuples[start].b;
        int end = start+1;

        while(end < R1Tuples.size() && R1Tuples[end].b == currB){
            end++;
        }

        std::vector<R1Tuple> bGroup(R1Tuples.begin() + start, R1Tuples.begin() + end);
        Ai.emplace_back(currB, bGroup);

        start = end;
    }

    std::vector<std::pair<int,std::vector<R2Tuple>>> Ci;
    start = 0;

    while(start < R2Tuples.size()){
        int currB = R2Tuples[start].b;
        int end = start+1;

        while(end < R2Tuples.size() && R2Tuples[end].b == currB){
            end++;
        }

        std::vector<R2Tuple> bGroup(R2Tuples.begin() + start, R2Tuples.begin() + end);
        Ci.emplace_back(currB, bGroup);

        start = end;
    }

    std::vector<ACpair> S;
    std::vector<ACpair> F;
    S.reserve(K_VAL);
    F.reserve(K_VAL);

    int i = 0, j = 0;
    while (i < (int)Ai.size() && j < (int)Ci.size()) {
        int b1 = Ai[i].first;
        int b2 = Ci[j].first;
        if (b1 == b2) {
            pointerSweep(Ai[i].second, Ci[j].second, p_val, K_VAL, S, F);
            i++;
            j++;
        } else if (b1 < b2) {
            i++;
        } else {
            j++;
        }
    }

    combine(S,F);

    // If S.size() == k, final estimate is k / p
    if ((int)S.size() == K_VAL) {
        double estimate = (double)K_VAL/p_val;
        std::cout << "Estimated join size: " << estimate << "\n";
    } else {
        // If we never filled S to k, the algorithm says z <= k^2
        std::cout << "Join size <= " << (K_VAL * (long long)K_VAL) << "\n";
    }

    return 0;

}