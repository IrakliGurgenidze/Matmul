#include <vector>

double EPSILON = 0.1;
int K_VAL = (int) (9.0 / (EPSILON*EPSILON));
double P_VAL = 1.0;

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

//TODO: pairwise hashing
double hash(int x) {
  return x;
}

double hashAC(double h1a, double h2c){
    return h1a-h2c;
}

int main() {

    std::vector<std::vector<int> > R1 = {{1, 0, 0}, {1, 1, 0}, {0, 0, 1}};
    std::vector<std::vector<int> > R2 = {{0, 1, 0}, {1, 0, 1}, {0, 0, 1}};

    std::vector<R1Tuple> R1Tuples;
    std::vector<R2Tuple> R2Tuples;

    // R1 -> (a,b,h1(a))
    for (int a = 0; a < R1.size(); a++) {
        for (int b = 0; b < R1[a].size(); b++) {
          if(R1[a][b] == 1) {
            double hashVal = hash(a);
            R1Tuples.push_back({a, b, hashVal});
          }
        }
    }

    //R2 -> (b,c,h2(c))
    for (int b = 0; b < R2.size(); b++) {
        for (int c = 0; c < R2[b].size(); c++) {
          if(R2[b][c] == 1) {
            double hashVal = hash(c);
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


    return 0;

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

    P_VAL = thresh;
}