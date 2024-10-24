


#include <iostream>
#include <vector>


#include "Solver.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/dimacs.hpp"

using namespace tempo;
class BIBD {
public:
  BIBD(Options &opt, const int v, const int b, const int r, const int k,
       const int lambda);

  Solver<> S;
  

  std::vector<BooleanVar<>> X;
  std::vector<BooleanVar<>> transpX;
  std::vector<BooleanVar<>> rowProduct;
  int v;
  int b;
  int r;
  int k;
  int lambda;
    
//    SubscriberHandle handlerToken;

  void print() {
    for (auto i{0}; i < v; ++i) {
      for (auto j{0}; j < b; ++j) {
        auto x{X[i * b + j]};
        std::cout << (S.boolean.isTrue(x) ? "1"
                                          : (S.boolean.isFalse(x) ? "0" : "."));
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    for (auto i{0}; i < v * (v - 1) / 2; ++i) {
      for (auto j{0}; j < b; ++j) {
        auto x{rowProduct[i * b + j]};
        std::cout << (S.boolean.isTrue(x) ? "1"
                                          : (S.boolean.isFalse(x) ? "0" : "."));
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
};

BIBD::BIBD(Options &opt, const int v, const int b, const int r, const int k,
           const int lambda)
    : S(opt), v(v), b(b), r(r), k(k), lambda(lambda)
{
    
    for (auto i{0}; i < v * b; ++i) {
        X.push_back(S.newBoolean());
    }

  for (auto j{0}; j < b; ++j) {
    for (auto i{0}; i < v; ++i) {
      transpX.push_back(X[i * b + j]);
    }
  }

  for (auto i{0}; i < v; ++i) {
    for (auto z{i + 1}; z < v; ++z) {
      for (auto j{0}; j < b; ++j) {
        rowProduct.push_back(X[i * b + j] and X[z * b + j]);
          rowProduct.back().extract(S);
      }
    }
  }

  std::vector<Literal<int>> scope;

  // rows
  auto x{X.begin()};
  for (auto i{0}; i < v; ++i) {
    S.postCardinality(x, x + b, true, r);
    S.postCardinality(x, x + b, false, b - r);
    x += b;
  }

  // columns
  x = transpX.begin();
  for (auto i{0}; i < b; ++i) {
    S.postCardinality(x, x + v, true, k);
    S.postCardinality(x, x + v, false, v - k);
    x += v;
  }

  auto y{rowProduct.begin()};
  // products
  for (auto i{0}; i < v * (v - 1) / 2; ++i) {
    S.postCardinality(y, y + b, true, lambda);
    S.postCardinality(y, y + b, false, b - lambda);
    y += b;
  }

  for (auto x : X)
    S.addToSearch(x);

  auto sat{S.satisfiable()};

  for (auto i{0}; i < v; ++i) {
    for (auto j{0}; j < b; ++j) {
      std::cout << S.boolean.value(X[i * b + j]);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for (auto i{0}; i < v * (v - 1) / 2; ++i) {
    for (auto j{0}; j < b; ++j) {
      std::cout << S.boolean.value(rowProduct[i * b + j]);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "result = " << sat << " #fails = " << S.num_fails << std::endl;
}

void PB() {
  Solver<float> solver;

  std::vector<BooleanVar<float>> X;
  for (auto i{0}; i < 10; ++i)
    X.push_back(solver.newBoolean());

  std::vector<Literal<float>> S{X[0] == true, X[1] == true, X[2] == false};
  std::vector<float> W{3.5, -5.2, 12.1};
  solver.postPseudoBoolean(S.begin(), S.end(), W.begin(), 7.1);

  S = {X[0] == true, X[3] == true, X[5] == false, X[7] == false};
  W = {1.5, 2.2, 2.1, 5.7};
  solver.postPseudoBoolean(S.begin(), S.end(), W.begin(), 7);

  S = {X[2] == true, X[4] == true, X[6] == false, X[8] == false, X[9] == false};
  W = {-1.5, -2.2, -2.1, -5.7, -3};
  solver.postPseudoBoolean(S.begin(), S.end(), W.begin(), -5);

  for (auto x : X)
    solver.addToSearch(x);

  auto sat{solver.satisfiable()};

  std::cout << "result = " << sat << " #fails = " << solver.num_fails
            << std::endl;

  for (auto x : X)
    std::cout << x << " " << solver.boolean.value(x) << std::endl;
}

int main() {

//    Options opt = tempo::parse(argc, argv);
    
    Options opt;
  //
      BIBD t1(opt, 7,7,3,3,1);
  //
      BIBD t2(opt, 7, 14, 6, 3, 2);
  //
    BIBD t3(opt, 6, 10, 5, 3, 2);
  //
  //      BIBD t4(opt, 46,69,9,6,1);

//  PB();
}
