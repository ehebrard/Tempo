/************************************************
 * Tempodisjunctive scheduling solver
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#include <iostream>
#include <vector>

#include "Solver.hpp"

using namespace tempo;

int main() {

  Solver<> S;

    auto makespan{S.newNumeric(0, 50)};
    auto origin{S.newConstant(0)};
    auto schedule{S.between(origin, makespan)};
//  auto schedule = S.newInterval(0, 50, 0, 0, 0, 50);

  auto s1{S.newNumeric()};
  auto t1{S.between(s1, s1 + 20)}; // x 5 = 100
  auto d1{S.newNumeric(5, 5)};

  auto s2{S.newNumeric()};
  auto t2{S.between(s2, s2 + 15)}; // x 10 = 150
  auto d2{S.newNumeric(10, 10)};

  auto s3{S.newNumeric()};
  auto t3{S.between(s3, s3 + 10)}; // x 12 = 120
  auto d3{S.newNumeric(12, 12)};

  S.post(s1 >= schedule.start);
  S.post(s2 >= 3);
  S.post(s3 >= 3);

  S.post(t1.end <= 25);
  S.post(t2.end <= 30);
  S.post(t3.end <= schedule.end);

  S.post(
      Cumulative<int>(schedule, S.newConstant(15), {t1, t2, t3}, {d1, d2, d3}));

  S.initializeSearch();
  S.propagate();

  std::cout << t1 << " (" << d1.min(S) << "): " << t1.start.min(S) << ".."
            << t1.end.max(S) << std::endl;

  std::cout << t2 << " (" << d2.min(S) << "): " << t2.start.min(S) << ".."
            << t2.end.max(S) << std::endl;

  std::cout << t3 << " (" << d3.min(S) << "): " << t3.start.min(S) << ".."
            << t3.end.max(S) << std::endl;
}
