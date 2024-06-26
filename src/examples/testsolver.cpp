/*************************************************************************
minicsp

Copyright 2010--2011 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/


#include <iostream>
#include <vector>


#include "Solver.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/dimacs.hpp"

using namespace tempo;

void test1(Options& opt) {

  Solver<float> S(opt);

  std::cout << S << std::endl;

  auto b0{S.newBoolean()};
  auto b1{S.newBoolean()};
  auto b2{S.newBoolean()};

  auto x0{S.newTemporal()};
  auto x1{S.newNumeric()};
  auto x2{S.newTemporal()};
  auto x3{S.newTemporal()};

  //    auto d0{S.newDisjunct(x0.before(x2,100), x0.after(x2,3))};
  auto d0{S.newDisjunct(x0.before(x2, 100), x2.before(x0, 3))};

  std::cout << S << std::endl;
  S.propagate();
  auto s1{S.saveState()};

  S.set(x0 <= float(18.999));

  S.set(x0 > float(-100));

  S.set(b1 == false);

  S.set(x1 < float(1000));

  S.set(x1 >= float(100));

  S.set(b2 == true);

  S.set(x3 > float(-9999999.99));

  std::cout << S << std::endl;
  S.propagate();
  auto s2{S.saveState()};

  S.set(b0 == true);

  S.set(x1 <= float(300));

  S.set(x0 >= float(0));

  std::cout << S << std::endl;
  S.propagate();
  auto s3{S.saveState()};

  S.set(d0 == true);

  std::cout << S << std::endl;

  std::cout << "restore to state " << s3 << std::endl;
  S.restoreState(s3);
  std::cout << S << std::endl;

  S.propagate();
  auto s4{S.saveState()};

  S.set(d0 == false);

  std::cout << S << std::endl;

  std::cout << "restore to state " << s4 << std::endl;
  S.restoreState(s4);
  std::cout << S << std::endl;

  std::cout << "restore to state " << s2 << std::endl;
  S.restoreState(s2);
  std::cout << S << std::endl;

  std::cout << "restore to state " << s1 << std::endl;
  S.restoreState(s1);
  std::cout << S << std::endl;
}

void test2(Options &opt) {

  Solver<> S(opt);

  auto schedule{S.newInterval(0, 100)};

  S.set(schedule.start.after(0));
  S.set(schedule.start.before(0));

  Interval<int> j0{S.newInterval(15, 15)};
  auto j1{S.newInterval(7, 10)};
  auto j2{S.newInterval(12, 12)};
  auto j3{S.newInterval(1, 5)};
  auto j4{S.newInterval(10, 10)};
  auto j5{S.newInterval(6, 6)};
  auto j6{S.newInterval(10, 12)};
  auto j7{S.newInterval(7, 7)};
  auto j8{S.newInterval(3, 3)};

  S.set(j0.start.after(schedule.start));
  S.set(j0.end.before(j1.start));
  S.set(j1.end.before(j2.start, 10));
  S.set(j2.end.before(schedule.end));

  S.set(j3.start.after(schedule.start));
  S.set(j3.end.before(j4.start));
  S.set(j4.end.before(j5.start));
  S.set(j5.end.before(schedule.end));

  S.set(j6.start.after(schedule.start));
  S.set(j6.end.before(j7.start));
  S.set(j7.end.before(j8.start));
  S.set(j8.end.before(schedule.end));

  DisjunctiveResource<int> R1({j0, j3, j6});
  DisjunctiveResource<int> R2({j1, j4, j7});
  DisjunctiveResource<int> R3({j2, j5, j8});

  std::vector<DisjunctVar<int>> X;

  R1.createOrderVariables(S, X);
  R2.createOrderVariables(S, X);
  R3.createOrderVariables(S, X);

  for (auto x : X)
    S.addToSearch(x);

  std::cout << S << std::endl;

  S.propagate();
  auto s1{S.saveState()};
  S.boolean_search_vars.remove_back(X[0]);
  S.set(X[0] == true);

  std::cout << S << std::endl;

  S.propagate();
  auto s2{S.saveState()};
  S.boolean_search_vars.remove_back(X[1]);
  S.set(X[1] == false);

  std::cout << S << std::endl;
  int s3;

  try {
    S.propagate();
    s3 = S.saveState();
    S.boolean_search_vars.remove_back(X[2]);
    S.set(X[2] == true);

    std::cout << S << std::endl;
  } catch (Failure<int> &f) {
    std::cout << "fail -> backtrack";
    S.restoreState(s3);

    std::cout << S << std::endl;

    s3 = S.saveState();
    S.boolean_search_vars.remove_back(X[2]);
    S.set(X[2] == false);

    std::cout << S << std::endl;

    S.restoreState(s2);

    std::cout << S << std::endl;

    S.restoreState(s1);

    std::cout << S << std::endl;
  }
}

void test3(Options &options) {

  Solver<> S(options);

  auto schedule{S.newInterval(0, 100)};

  S.set(schedule.start.after(0));
  S.set(schedule.start.before(0));

  Interval<int> j0{S.newInterval(15, 15)};
  auto j1{S.newInterval(7, 10)};
  auto j2{S.newInterval(12, 12)};
  auto j3{S.newInterval(1, 5)};
  auto j4{S.newInterval(10, 10)};
  auto j5{S.newInterval(6, 6)};
  auto j6{S.newInterval(10, 12)};
  auto j7{S.newInterval(7, 7)};
  auto j8{S.newInterval(3, 3)};

  S.set(j0.start.after(schedule.start));
  S.set(j0.end.before(j1.start));
  S.set(j1.end.before(j2.start, 10));
  S.set(j2.end.before(schedule.end));

  S.set(j3.start.after(schedule.start));
  S.set(j3.end.before(j4.start));
  S.set(j4.end.before(j5.start));
  S.set(j5.end.before(schedule.end));

  S.set(j6.start.after(schedule.start));
  S.set(j6.end.before(j7.start));
  S.set(j7.end.before(j8.start));
  S.set(j8.end.before(schedule.end));

  DisjunctiveResource<int> R1({j0, j3, j6});
  DisjunctiveResource<int> R2({j1, j4, j7});
  DisjunctiveResource<int> R3({j2, j5, j8});

  std::vector<DisjunctVar<int>> X;

  R1.createOrderVariables(S, X);
  R2.createOrderVariables(S, X);
  R3.createOrderVariables(S, X);

  for (auto x : X)
    S.addToSearch(x);

  //  std::cout << S << std::endl;

  MinimizationObjective<int, TemporalVar<int>> makespan(schedule.end);

  S.optimize(makespan);

  //  std::cout << duration.value() << std::endl;
}

void test4(Options &opt) {

  Solver<> S(opt);

  std::vector<BooleanVar<>> X;

  dimacs::parse(opt.instance_file, S, X);

  for (auto x : X)
    S.addToSearch(x);

  //  std::cout << S.clauses << std::endl;

  auto sat{S.satisfiable()};

  std::cout << "result = " << sat << " #fails = " << S.num_fails << std::endl;
}

void test5(Options &opt) {

  Solver<> S(opt);

  auto schedule{S.newInterval()};
  std::vector<DisjunctiveResource<>> resources;
  std::vector<Interval<>> Intervals;

  osp::parse(opt.instance_file, S, schedule, Intervals, resources);

  std::vector<BooleanVar<>> X;
  for (auto &R : resources) {
    //        auto disjuncts{X.end()};
    auto s{X.size()};
    R.createOrderVariables(S, X);
    if (opt.edge_finding) {
      S.postEdgeFinding(schedule, R.begin(), R.end(), X.begin() + s);
    }
    if (opt.transitivity) {
      S.postTransitivity(schedule, R.begin(), R.end(), X.begin() + s);
    }
  }

  for (auto x : X)
    S.addToSearch(x);

  //    std::cout << S << std::endl;

  MinimizationObjective<int, TemporalVar<int>> makespan(schedule.end);

  S.set(schedule.end.before(opt.ub));

  S.optimize(makespan);
}

class BIBD {
public:
  BIBD(Options &opt, const int v, const int b, const int r, const int k,
       const int lambda);

  Solver<> S;
  

  std::vector<BooleanVar<>> X;
  std::vector<BooleanVar<>> transpX;
  std::vector<BooleanExpression<>> rowProduct;
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
//      ,handlerToken(S.ChoicePoint.subscribe_handled([this]() { this->print(); }))
{

  //  Solver<> S(opt);

  for (auto i{0}; i < v * b; ++i)
    X.push_back(S.newBoolean());

  for (auto j{0}; j < b; ++j) {
    for (auto i{0}; i < v; ++i) {
      transpX.push_back(X[i * b + j]);
    }
  }

  for (auto i{0}; i < v; ++i) {
    for (auto z{i + 1}; z < v; ++z) {
      for (auto j{0}; j < b; ++j) {
        rowProduct.push_back(X[i * b + j] && X[z * b + j]);
        rowProduct.back().extract(S);
      }
    }
  }

  //    std::cout << "rows\n";

  // rows
  auto x{X.begin()};
  for (auto i{0}; i < v; ++i) {
    S.postCardinality(x, x + b, true, r);
    S.postCardinality(x, x + b, false, b - r);
    x += b;
  }

  //    std::cout << "columns\n";

  // columns
  x = transpX.begin();
  for (auto i{0}; i < b; ++i) {
    S.postCardinality(x, x + v, true, k);
    S.postCardinality(x, x + v, false, v - k);
    x += v;
  }

  //    std::cout << "products\n";

  auto y{rowProduct.begin()};
  // products
  for (auto i{0}; i < v * (v - 1) / 2; ++i) {

    //        std::cout << "product card " << i << std::endl;

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

int main(int argc, char *argv[]) {

  Options opt = tempo::parse(argc, argv);

    BIBD t1(opt, 7,7,3,3,1);

    BIBD t2(opt, 7, 14, 6, 3, 2);

  BIBD t3(opt, 6, 10, 5, 3, 2);

      BIBD t4(opt, 46,69,9,6,1);
}
