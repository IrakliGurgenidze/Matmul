#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>


double EPSILON = 0.1;
int K_VAL = (int) (9.0 / (EPSILON*EPSILON));

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

double p_val = 1.0;

const uint64_t PRIME = 4294967311ULL;

uint64_t a1,b1,a2,b2;

void initPairwiseHashes(){
    std::random_device rd;
    std::mt19937_64 gen((rd()));
    std::uniform_int_distribution<uint64_t> dis(1, PRIME-1);
    a1 = dis(gen);
    b1 = dis(gen);
    a2 = dis(gen);
    b2 = dis(gen);
}

//TODO: pairwise hashing
double hash(int x, uint64_t a, uint64_t b) {
  uint64_t hash_val = (a * static_cast<uint64_t>(x) + b) % PRIME;
  return static_cast<double> (hash_val) / static_cast<double> (PRIME);
}

double hashAC(double h1a, double h2c){
    return (h1a-h2c) - floor(h1a-h2c);
}

//TODO: Make combine return similar to psuedocode?
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

int main() {

    int n = 100;
    double density = 0.1; // probability that an entry is 1

    // Create two 100x100 matrices initialized with 0.
    std::vector<std::vector<int>> R1(n, std::vector<int>(n, 0));
    std::vector<std::vector<int>> R2(n, std::vector<int>(n, 0));

    // Set up a random generator:
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(density);  // each entry is 1 with probability 'density'

    // Fill R1 randomly
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            R1[i][j] = d(gen) ? 1 : 0;
        }
    }

    // Fill R2 randomly
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            R2[i][j] = d(gen) ? 1 : 0;
        }
    }


    std::vector<R1Tuple> R1Tuples;
    std::vector<R2Tuple> R2Tuples;

    initPairwiseHashes();

    // R1 -> (a,b,h1(a))
    for (int a = 0; a < R1.size(); a++) {
        for (int b = 0; b < R1[a].size(); b++) {
          if(R1[a][b] == 1) {
            double hashVal = hash(a, a1, b1);
            R1Tuples.push_back({a, b, hashVal});
          }
        }
    }

    //R2 -> (b,c,h2(c))
    for (int b = 0; b < R2.size(); b++) {
        for (int c = 0; c < R2[b].size(); c++) {
          if(R2[b][c] == 1) {
            double hashVal = hash(c, a2, b2);
            R2Tuples.push_back({b, c, hashVal});
          }
        }
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
        Ai.push_back({currB, bGroup});

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
        Ci.push_back({currB, bGroup});

        start = end;
    }


    //TODO: Start psuedocode process to estimate size now that we have Ai and Ci sorted properly

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