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
#include <sstream>
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
template <typename T> class Interval;

struct PrecedenceEncoding {
  int x;
  int y;
  bool strict;
};

template <typename T> class Incrementality : public Constraint<T> {

private:
  Solver<T> &solver;
  std::vector<Interval<T>> task;

  // precedence[i][j] <=> (e_i <= s_j or e_i > s_j)
  Matrix<Literal<T>> precedence;

public:
  // tasks that are ordered with respect to all other tasks AND are only after
  // tasks in front are in front tasks that are ordered with respect to all
  // other tasks AND are only before tasks in back are in back other tasks are
  // in free_tasks

  SparseSet<int, Reversible<size_t>> free_tasks;
  Reversible<T> lb;
  Reversible<T> ub;

private:
  //    int lb_witness{-1};
  //    int ub_witness{-1};

  std::vector<Reversible<size_t>> numRanked;

  std::vector<PrecedenceEncoding> trigger;

public:
  template <concepts::typed_range<Interval<T>> Tasks>
  Incrementality(Solver<T> &solver, Tasks &&tasks, Matrix<Literal<T>> precs);
  virtual ~Incrementality();

  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  std::string asciiArt(const int i) const;

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DBG_INCRP
  int debug_flag{2};
  size_t count_ordered(const int i);
  void verify(const char *msg);
#endif
};

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks>
Incrementality<T>::Incrementality(Solver<T> &solver, Tasks &&tasks,
                                  Matrix<Literal<T>> precs)
    : solver(solver), task(std::forward<Tasks>(tasks).begin(),
                           std::forward<Tasks>(tasks).end()),
      precedence(std::move(precs)),
      free_tasks(precedence.numRows(), &(solver.getEnv())),
      lb(-Constant::Infinity<T>, &(solver.getEnv())),
      ub(Constant::Infinity<T>, &(solver.getEnv())) {

  Constraint<T>::priority = Priority::High;

  free_tasks.fill();

  T minest{Constant::Infinity<T>};
  T maxlct{-Constant::Infinity<T>};

  int n{static_cast<int>(precedence.numRows())};
  for (int i{0}; i < n; ++i) {
    numRanked.emplace_back(0, &solver.getEnv());
    //      if(est(i) < minest) {
    //          minest = est(i);
    //          lb_witness = i;
    //      }
    //      if(lct(i) > maxlct) {
    //          maxlct = lct(i);
    //          ub_witness = i;
    //      }
    minest = std::min(minest, est(i));
    maxlct = std::max(maxlct, lct(i));
  }
  lb = minest;
  ub = maxlct;
}

template <typename T> Incrementality<T>::~Incrementality() {}

template <typename T> T Incrementality<T>::est(const unsigned i) const {
  return task[i].getEarliestStart(solver);
}

template <typename T> T Incrementality<T>::lst(const unsigned i) const {
  return task[i].getLatestStart(solver);
}

template <typename T> T Incrementality<T>::ect(const unsigned i) const {
  return task[i].getEarliestEnd(solver);
}

template <typename T> T Incrementality<T>::lct(const unsigned i) const {
  return task[i].getLatestEnd(solver);
}

template <typename T> T Incrementality<T>::minduration(const unsigned i) const {
  return task[i].minDuration(solver);
}

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
        if (precedence(i, j).variable() == precedence(j, i).variable()) {
          if (i < j) {
            solver.wake_me_on(precedence(i, j), this->id());
            trigger.emplace_back(i, j, true);
            solver.wake_me_on(precedence(j, i), this->id());
            trigger.emplace_back(j, i, true);
          }
        } else {
          solver.wake_me_on(precedence(i, j), this->id());
          trigger.emplace_back(i, j, false);
          solver.wake_me_on(~precedence(i, j), this->id());
          trigger.emplace_back(j, i, true);
        }
      }
    }
  }
}

template <typename T>
bool Incrementality<T>::notify(const Literal<T>
#ifdef DBG_INCR
                                   l
#endif
                               ,
                               const int r) {

  //    auto n{precedence.size()};
  //    auto i{r / n};
  //    auto j{r % n};

  auto i{trigger[r].x};
  auto j{trigger[r].y};
  auto strict{trigger[r].strict};

#ifdef DBG_INCR
  if (DBG_INCR) {
    //    std::cout << "\ntrigger " << r << std::endl;
    std::cout << "\nt" << std::left << std::setw(3) << i << ": " << asciiArt(i)
              << std::endl;
    std::cout << "t" << std::left << std::setw(3) << j << ": " << asciiArt(j)
              << std::endl;

    std::cout << "[" << this->id() << "] notify " << solver.pretty(l)
              << ", i.e., " << (i)
              << (strict ? " is strictly after " : " is not before ") << (j);
  }
#endif

  //    if(precedence(i,j).variable() == precedence(j,i).variable())
  if (strict or (solver.boolean.satisfied(precedence(j, i)) and
                 (solver.propagationStamp(precedence(i, j)) >
                  solver.propagationStamp(precedence(j, i))))) {

#ifdef DBG_INCR
    if (DBG_INCR) {
      std::cout << " ==> ordered" << (strict ? "" : "*") << "!\n";
    }
#endif

    ++numRanked[i];
    ++numRanked[j];
  }
#ifdef DBG_INCR
  else if (DBG_INCR) {
    std::cout << " ==> not completely defined\n";
  }
#endif

  auto n{precedence.numRows() - 1};

#ifdef DBG_INCR
  if (DBG_INCR) {
    auto ci{count_ordered(i)};
    auto cj{count_ordered(j)};

    std::cout << numRanked[i] << "|" << ci << " & " << numRanked[j] << "|" << cj
              << " / " << n << " " << free_tasks.has(i) << free_tasks.has(j)
              << std::endl;
  }
#endif

  bool rmtask{false};

  if (free_tasks.has(i) and numRanked[i] == n) {
    free_tasks.remove_back(i);

    //      if(i)

#ifdef DBG_INCR
    if (DBG_INCR) {
      std::cout << "rm " << i << " from free tasks :" << free_tasks
                << std::endl;
    }
#endif

    rmtask = true;
  }

  if (free_tasks.has(j) and numRanked[j] == n) {
    free_tasks.remove_back(j);

#ifdef DBG_INCR
    if (DBG_INCR) {
      std::cout << "rm " << j << " from free tasks :" << free_tasks
                << std::endl;
    }
#endif

    rmtask = true;
  }

  //    std::cout << "return " << rmtask << std::endl;

  return rmtask;

  //  return false;
}

#ifdef DBG_INCRP

template <typename T> size_t Incrementality<T>::count_ordered(const int i) {
  auto n{static_cast<int>(precedence.numRows()) - 1};
  size_t num_ordered{0};
  for (int j{0}; j <= n; ++j) {
    if (i != j) {
      auto ordered{((solver.boolean.falsified(precedence(i, j)) or
                     solver.boolean.satisfied(precedence(i, j))) and
                    (solver.boolean.falsified(precedence(j, i)) or
                     solver.boolean.satisfied(precedence(j, i))))};
      if (ordered) {
        ++num_ordered;
      }
    }
  }
  return num_ordered;
}

template <typename T> void Incrementality<T>::verify(const char *msg) {
  auto n{precedence.numRows() - 1};
  for (unsigned i{0}; i <= n; ++i) {
    //    std::cout << std::setw(2) << i << ": " << std::setw(3)
    //              << (n - numRanked[i]);
    //        unsigned num_ordered{0};
    //        for (unsigned j{0}; j <= n; ++j) {
    //
    //            //      if (i == j) {
    //            //        std::cout << " *";
    //            //      } else
    //
    //            if(i != j) {
    //                auto ordered{((solver.boolean.falsified(precedence(i, j))
    //                or
    //                               solver.boolean.satisfied(precedence(i, j)))
    //                               and
    //                              (solver.boolean.falsified(precedence(j, i))
    //                              or
    //                               solver.boolean.satisfied(precedence(j,
    //                               i))))};
    //                if (ordered) {
    //                    ++num_ordered;
    //                    //          std::cout << " " << i << "." << j;
    //
    //                    //                    if((i==27 and j==29) or ()) {
    //
    //                }
    //
    //
    //
    //                //                  num_ordered += ordered;
    //            }
    //            //                                std::cout << (i==j or
    //            ordered);
    //        }
    size_t num_ordered{count_ordered(i)};
    //    std::cout << std::endl;

    if (num_ordered != numRanked[i]) {
      std::cout << msg << " [" << this->id() << "] discrepancy: " << num_ordered
                << "/" << numRanked[i] << std::endl;

      for (unsigned j{0}; j <= n; ++j) {

        if (i == j) {
          std::cout << " *\n";
        } else {
          auto ordered{((solver.boolean.falsified(precedence(i, j)) or
                         solver.boolean.satisfied(precedence(i, j))) and
                        (solver.boolean.falsified(precedence(j, i)) or
                         solver.boolean.satisfied(precedence(j, i))))};

          std::cout << "\n " << i << "." << j << " -->\n";
          if (ordered) {
            if (solver.boolean.falsified(precedence(i, j))) {
              std::cout << solver.pretty(precedence(i, j)) << " is falsified\n";
            }
            if (solver.boolean.satisfied(precedence(i, j))) {
              std::cout << solver.pretty(precedence(i, j)) << " is satisfied\n";
            }

            if (solver.boolean.falsified(precedence(j, i))) {
              std::cout << solver.pretty(precedence(j, i)) << " is falsified\n";
            }
            if (solver.boolean.satisfied(precedence(j, i))) {
              std::cout << solver.pretty(precedence(j, i)) << " is satisfied\n";
            }
          } else {
            if (not solver.boolean.falsified(precedence(i, j))) {
              std::cout << solver.pretty(precedence(i, j))
                        << " is not falsified\n";
            }
            if (not solver.boolean.satisfied(precedence(i, j))) {
              std::cout << solver.pretty(precedence(i, j))
                        << " is not satisfied\n";
            }

            if (not solver.boolean.falsified(precedence(j, i))) {
              std::cout << solver.pretty(precedence(j, i))
                        << " is not falsified\n";
            }
            if (not solver.boolean.satisfied(precedence(j, i))) {
              std::cout << solver.pretty(precedence(j, i))
                        << " is not satisfied\n";
            }
          }
          std::cout << "[" << est(i) << "-" << ect(i) << ".." << lst(i) << "-"
                    << lct(i) << "] vs [" << est(j) << "-" << ect(j) << ".."
                    << lst(j) << "-" << lct(j) << "]\n";
        }
      }
      std::cout << std::endl;

      exit(1);
    }
  }
}
#endif

template <typename T> void Incrementality<T>::propagate() {

  T minest{Constant::Infinity<T>};
  T maxlct{-Constant::Infinity<T>};

  //    auto dbg{false;}
  //    if(free_tasks.size() <= free_tasks.capacity() / 2) {
  //        dbg = true;
  //
  //    }

  for (auto i : free_tasks) {
    minest = std::min(minest, est(i));
    maxlct = std::max(maxlct, lct(i));
    //
    //        if(dbg) {
    //            std::cout << "[" << est(i) << ".." << lct(i) << "]"
    //        }
  }

  if (lb < minest)
    lb = minest;

  if (ub > maxlct)
    ub = maxlct;

#ifdef DBG_INCRP
  std::stringstream ss;
  ss << "propagate @lvl" << solver.level() << std::endl;
  verify(ss.str().c_str());
#endif
}

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

template <typename T>
std::string Incrementality<T>::asciiArt(const int i) const {
  std::stringstream ss;
  ss
      //    << std::setw(3) << std::right << mindemand(i) << "x"
      << std::setw(3) << std::left << minduration(i) << " " << std::right;
  for (auto k{0}; k < est(i); ++k) {
    ss << " ";
  }
  auto est_i{est(i)};
  auto ect_i{ect(i)};

  if (est_i == -Constant::Infinity<T>) {
    ss << "...";
    est_i = -1;
    ect_i = minduration(i);
  } else {
    ss << "[";
  }
  for (auto k{est_i + 1}; k < ect(i); ++k) {
    ss << "=";
  }
  if (ect_i < lct(i))
    ss << "|";
  if (lct(i) == Constant::Infinity<T>) {
    ss << "... " << est(i) << "...";
  } else {
    for (auto k{ect_i + 1}; k < lct(i); ++k) {
      ss << ".";
    }
    ss << "] " << est(i) << "-" << ect(i) << ".." << lct(i);
  }
  return ss.str();
}

} // namespace tempo

#endif
