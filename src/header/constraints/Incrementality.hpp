/************************************************
 * Tempo Incrementality.hpp
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

#ifndef TEMPO_INCREMENTALITY_HPP
#define TEMPO_INCREMENTALITY_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Model.hpp"
#include "ReversibleObject.hpp"
#include "constraints/Constraint.hpp"
#include "util/Matrix.hpp"
#include "util/SparseSet.hpp"
#include "util/traits.hpp"

namespace tempo {

template <typename T> class Solver;

template <typename T> class Incrementality : public Constraint<T> {
private:
  Solver<T> &solver;

  // precedence[i][j] <=> (e_i <= s_j or e_i > s_j)
  Matrix<Literal<T>> precedence;

  // tasks that are ordered with respect to all other tasks AND are only after
  // tasks in front are in front tasks that are ordered with respect to all
  // other tasks AND are only before tasks in back are in back other tasks are
  // in free_tasks
  SparseSet<int, Reversible<size_t>> free_tasks;

  std::vector<Reversible<size_t>> numRanked;

  std::vector<std::pair<int, int>> scope;

public:
  Incrementality(Solver<T> &solver, Matrix<Literal<T>> precs);
  virtual ~Incrementality();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
Incrementality<T>::Incrementality(Solver<T> &solver, Matrix<Literal<T>> precs)
    : solver(solver), precedence(std::move(precs)),
      free_tasks(precs.numRows(), &(solver.getEnv())) {

  Constraint<T>::priority = Priority::High;

  //              for(auto : precedence) {
  //                  numRanked.emplace_back(0, &solver);
  //              }

  for (size_t i{0}; i < precedence.numRows(); ++i) {
    numRanked.emplace_back(0, &solver.getEnv());
  }
}

template <typename T> Incrementality<T>::~Incrementality() {}

template <typename T> void Incrementality<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < precedence.numRows(); ++i) {
    for (size_t j{0}; j < precedence.numColumns(); ++j) {
      if (i != j) {
        solver.wake_me_on(precedence(i, j), this->id());
        scope.emplace_back(i, j);

        std::cout << " " << solver.pretty(precedence(i, j));
      } else
        std::cout << " *** ";
    }
    std::cout << std::endl;
  }
}

template <typename T>
bool Incrementality<T>::notify(const Literal<T> l, const int r) {

  //    auto n{precedence.size()};
  //    auto i{r / n};
  //    auto j{r % n};

  auto i{scope[r].first};
  auto j{scope[r].second};

  std::cout << "notify " << solver.pretty(l) << " / check "
            << solver.pretty(precedence(j, i))
            << (solver.boolean.satisfied(precedence(j, i))
                    ? " t"
                    : (solver.boolean.falsified(precedence(j, i)) ? " f"
                                                                  : " u"))
            << std::endl;

  //    if(l.sign())
  std::cout << (j) << " is not before " << (i);
  //    if(precedence(i,j).variable() == precedence(j,i).variable())
  if (solver.boolean.satisfied(precedence(j, i))) {
    std::cout << " ==> overlap\n";
    ++numRanked[i];
    ++numRanked[j];

  } else if (solver.boolean.falsified(precedence(j, i))) {
    std::cout << " ==> precedence\n";
    ++numRanked[i];
    ++numRanked[j];
  }

  auto n{precedence.numRows() - 1};

  std::cout << numRanked[i] << " & " << numRanked[j] << " / "
            << precedence.numRows() - 1 << std::endl;

  return false;
}

template <typename T> void Incrementality<T>::propagate() {}

template <typename T>
void Incrementality<T>::xplain(const Literal<T>, const hint,
                               std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &Incrementality<T>::display(std::ostream &os) const {
  os << "Incrementality";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif

  return os;
}

template <typename T>
std::ostream &Incrementality<T>::print_reason(std::ostream &os,
                                              const hint) const {
  os << "Incrementality";
  return os;
}

} // namespace tempo

#endif
