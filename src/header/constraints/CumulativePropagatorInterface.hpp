/************************************************
 * Tempo CumulativePropagatorInterface.hpp
 * Interface for cumulative constraint propagators.
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

#ifndef TEMPO_CUMULATIVEINTERFACE_HPP
#define TEMPO_CUMULATIVEINTERFACE_HPP

#include <cassert>
#include <map>
#include <vector>
#include <sstream>
#include <numeric>


#include "Explanation.hpp"
#include "Global.hpp"
#include "Model.hpp"
#include "constraints/Constraint.hpp"
#include "util/List.hpp"


namespace tempo {



template<typename T>
class Solver;


template <typename T> class CumulativePropagatorInterface : public Constraint<T> {
protected:
  Solver<T> &solver;
  NumericVar<T> capacity;

  std::vector<Interval<T>> task;
    std::vector<NumericVar<T>> demand;
  bool sign{true};
    
    std::vector<Literal<T>> pruning;
    std::vector<std::vector<Literal<T>>> explanation;
    Reversible<size_t> num_explanations;
    

public:
  template <typename ItTask, typename ItNVar>
    CumulativePropagatorInterface(Solver<T> &solver, const NumericVar<T> capacity,
                        const ItTask beg_task, const ItTask end_task,
                        const ItNVar beg_dem);

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;
    T mindemand(const unsigned i) const;

  void post(const int idx) override;
    bool notify(const Literal<T>, const int rank) override;
    void doPruning();

  // create a new explanation for a new pruning
  hint newExplanation();
  // function used in explanation
  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;


  std::ostream &display(std::ostream &os) const override;
  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;
};


template <typename T>
std::string CumulativePropagatorInterface<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << "t" << task[i].id() << ": ["
     << (sign == bound::lower ? est(i) : -lct(i)) << ".."
     << (sign == bound::lower ? lct(i) : -est(i)) << "] (p=" << minduration(i)
     << ", c=" << mindemand(i) << ")";
  return ss.str();
}

template <typename T>
std::string CumulativePropagatorInterface<T>::asciiArt(const int i) const {
  std::stringstream ss;
  ss << std::setw(3) << std::right << mindemand(i) << "x" << std::setw(3)
     << std::left << minduration(i) << " " << std::right;
  for (auto k{0}; k < est(i); ++k) {
    ss << " ";
  }
  ss << "[";
  for (auto k{est(i) + 1}; k < ect(i); ++k) {
    ss << "=";
  }
  if (lct(i) == Constant::Infinity<T>) {
    ss << "=... " << est(i) << "...";
  } else {
    for (auto k{ect(i)}; k < lct(i); ++k) {
      ss << ".";
    }
    ss << "] " << est(i) << ".." << lct(i);
  }
  return ss.str();
}

template <typename T> T CumulativePropagatorInterface<T>::est(const unsigned i) const {
  return (sign == bound::lower ? task[i].getEarliestStart(solver)
                               : -task[i].getLatestEnd(solver));
}

template <typename T> T CumulativePropagatorInterface<T>::lst(const unsigned i) const {
  return (sign == bound::lower ? task[i].getLatestStart(solver)
                               : -task[i].getEarliestEnd(solver));
}

template <typename T> T CumulativePropagatorInterface<T>::ect(const unsigned i) const {
  return (sign == bound::lower ? task[i].getEarliestEnd(solver)
                               : -task[i].getLatestStart(solver));
}

template <typename T> T CumulativePropagatorInterface<T>::lct(const unsigned i) const {
  return (sign == bound::lower ? task[i].getLatestEnd(solver)
                               : -task[i].getEarliestStart(solver));
}

template <typename T>
T CumulativePropagatorInterface<T>::minduration(const unsigned i) const {
  return task[i].minDuration(solver);
}

template <typename T>
T CumulativePropagatorInterface<T>::maxduration(const unsigned i) const {
  return task[i].maxDuration(solver);
}

template <typename T>
T CumulativePropagatorInterface<T>::mindemand(const unsigned i) const {
  return demand[i].min(solver);
}

//template <typename T>
//T CumulativeTimetabling<T>::maxdemand(const unsigned i) const {
//  return demand[i].max(solver);
//}

//template <typename T>
//T CumulativeTimetabling<T>::energy(const unsigned i) const {
//  return mindemand(i) * minduration(i);
//}

template <typename T>
template <typename ItTask, typename ItNVar>
CumulativePropagatorInterface<T>::CumulativePropagatorInterface(Solver<T> &solver,
                                                const NumericVar<T> cap,
                                                const ItTask beg_task,
                                                const ItTask end_task,
                                                                  const ItNVar beg_dem)
    : solver(solver), num_explanations(0, &(solver.getEnv())) {
  capacity = cap;

//  Constraint<T>::priority = Priority::Medium;

        auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    task.push_back(*jp);
      demand.push_back(*dp);
      ++dp;
  }
}


template <typename T> void CumulativePropagatorInterface<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_CEDGEFINDING
  if (DBG_CEDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < task.size(); ++i) {
    solver.wake_me_on(lb<T>(task[i].start.id()), this->id());
    solver.wake_me_on(ub<T>(task[i].end.id()), this->id());
    solver.wake_me_on(lb<T>(demand[i].id()), this->id());
  }
}

template <typename T>
bool CumulativePropagatorInterface<T>::notify(const Literal<T>, const int) {
  return true;
}



template <typename T> hint CumulativePropagatorInterface<T>::newExplanation() {
  auto e_idx{num_explanations};
  if (explanation.size() <= e_idx) {
    explanation.resize(e_idx + 1);
  } else {
    explanation[e_idx].clear();
  }
  ++num_explanations;
  return static_cast<hint>(e_idx);
}



template <typename T> void CumulativePropagatorInterface<T>::doPruning() {
    
#ifdef STATS_TT
    num_pruning += pruning.size();
    if(not pruning.empty()) {
        ++num_useful;
        std::cout << "TT" << this->id() << ": " << static_cast<double>(num_useful)/static_cast<double>(num_prop) << " (" << num_pruning << ")\n";
    }
#endif
    

#ifdef DBG_TT
  if (DBG_TT) {
    if (not pruning.empty())
      std::cout << "apply pruning" << std::endl;
  }
#endif

  auto h{static_cast<hint>(num_explanations - pruning.size())};
  for (auto p : pruning) {
#ifdef DBG_TT
    if (DBG_TT) {
      std::cout << "pruning (" << this->id() << "/"
                << solver.num_cons_propagations << "): " << p << std::endl;
    }
#endif

    solver.set(p, {this, h++});
  }
  pruning.clear();
}

template <typename T>
void CumulativePropagatorInterface<T>::xplain(const Literal<T>, const hint h,
                                      std::vector<Literal<T>> &Cl) {

  for (auto p : explanation[h]) {
    Cl.push_back(p);
  }
}


template <typename T>
std::ostream &CumulativePropagatorInterface<T>::display(std::ostream &os) const {
  os << "Cumulative Time-Tabling";

#ifdef DBG_CEDGEFINDING
  os << "[" << this->id() << "]";
#endif

  os << "(";
  for (auto &t : task) {
    std::cout << " t" << t.id();
  }
  std::cout << " )";
  return os;
}

} // namespace tempo

#endif
