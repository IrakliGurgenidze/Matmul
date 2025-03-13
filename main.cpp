#include <stdlib>

struct R1Tuple {
  int a;
  int b;
  double h1a;
}

struct R2Tuple {
  int b;
  int c;
  double h2c;
}

//TODO: pairwise hashing
double hash(int x) {
  return x;
}

int main() {

  std::vector<vector<int>> R1 = {{1, 0, 0}, {1, 1, 0}, {0, 0, 1}};
  std::vector<vector<int>> R2 = {{0, 1, 0}, {1, 0, 1}, {0, 0, 1}};

  std::vector<R1Tuple> R1Tuples;
  std::vector<R2Tuple> R2Tuples;

  for (int a = 0; a < R1.size(); a++) {
    for (int b = 0; b < R1[a].size(); b++) {
      if(R1[a][b] == 1) {
        double hashVal = hash(a);
        R1Tuples.push_back({a, b, hashVal});
      }
    }
  }

  for (int b = 0; b < R1.size(); b++) {
    for (int c = 0; c < R1[b].size(); c++) {
      if(R1[b][c] == 1) {
        double hashVal = hash(c);
        R1Tuples.push_back({b, c, hashVal});
      }
    }
  }

  return 0;

}