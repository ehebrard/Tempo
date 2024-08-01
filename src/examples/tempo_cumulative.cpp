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

  auto schedule = S.newInterval(0, Constant::Infinity<int>, 0, 0);

  auto t1 = S.newInterval(10, 10, 0, 20);
  auto d1 = S.newConstant(3);

  auto t2 = S.newInterval(15, 15, 0, 15);
  auto d2 = S.newConstant(1);

  auto t3 = S.newInterval(7, 7, 0, 20);
  auto d3 = S.newConstant(2);

  auto t4 = S.newInterval(13, 13, 5, 17);
  auto d4 = S.newConstant(2);

  auto t5 = S.newInterval(12, 12, 2, 17);
  auto d5 = S.newConstant(1);

  auto t6 = S.newInterval(11, 11, 5, 19);
  auto d6 = S.newConstant(1);

  auto t7 = S.newInterval(5, 5, 10, 25);
  auto d7 = S.newConstant(2);

  //    std::cout << t1 << std::endl;
  //    std::cout << t2 << std::endl;
  //    std::cout << t3 << std::endl;
  //    std::cout << t4 << std::endl;
  //    std::cout << t5 << std::endl;
  //    std::cout << t6 << std::endl;
  //    std::cout << t7 << std::endl;
  //    std::cout << S << std::endl;

  S.post(t1.end <= schedule.end);
  S.post(t2.end <= schedule.end);
  S.post(t3.end <= schedule.end);
  S.post(t4.end <= schedule.end);
  S.post(t5.end <= schedule.end);
  S.post(t6.end <= schedule.end);
  S.post(t7.end <= schedule.end);

  S.post(Cumulative<int>(S.newConstant(4), {t1, t2, t3, t4, t5, t6, t7},
                         {d1, d2, d3, d4, d5, d6, d7}));

  //    std::cout << S << std::endl;

  //    auto sat{S.satisfiable()};

  S.minimize(schedule.duration);

  //    std::cout << sat << std::endl;

  std::cout << t1 << " (" << S.numeric.lower(d1)
            << "): " << S.numeric.lower(t1.start) << ".."
            << S.numeric.upper(t1.end) << std::endl;
  std::cout << t2 << " (" << S.numeric.lower(d2)
            << "): " << S.numeric.lower(t2.start) << ".."
            << S.numeric.upper(t2.end) << std::endl;
  std::cout << t3 << " (" << S.numeric.lower(d3)
            << "): " << S.numeric.lower(t3.start) << ".."
            << S.numeric.upper(t3.end) << std::endl;
  std::cout << t4 << " (" << S.numeric.lower(d4)
            << "): " << S.numeric.lower(t4.start) << ".."
            << S.numeric.upper(t4.end) << std::endl;
  std::cout << t5 << " (" << S.numeric.lower(d5)
            << "): " << S.numeric.lower(t5.start) << ".."
            << S.numeric.upper(t5.end) << std::endl;
  std::cout << t6 << " (" << S.numeric.lower(d6)
            << "): " << S.numeric.lower(t6.start) << ".."
            << S.numeric.upper(t6.end) << std::endl;
  std::cout << t7 << " (" << S.numeric.lower(d7)
            << "): " << S.numeric.lower(t7.start) << ".."
            << S.numeric.upper(t7.end) << std::endl;
}
