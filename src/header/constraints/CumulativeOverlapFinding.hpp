/************************************************
 * Tempo CumulativeOverlapFinding.hpp
 * Implementation of the "strong edge-finding" algorithm as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder
 *rule for the cumulative constraint. In Proceedings of the Twenty-Fifth
 *International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI
 *Press, 3103–3109.
 *
 * Copyright 2024 Emmanuel Hebrard, Roger Kameugne
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
//
// Created by Roger Kameugne on 12/11/2024.
//



#ifndef CUMULATIVEOVERLAPFINDING_HPP
#define CUMULATIVEOVERLAPFINDING_HPP

#include <cassert>
#include <map>
#include <vector>
#include <sstream>
#include <numeric>


#include "Explanation.hpp"
#include "Global.hpp"
#include "Model.hpp"
#include "constraints/CumulativeEdgeFinding.hpp"
#include "util/List.hpp"
namespace tempo {

//template <typename T = int> struct Timepoint {
//
//  Timepoint() {}
//  Timepoint(T time, T increment, T incrementMax)
//      : time(time), increment(increment), incrementMax(incrementMax) {}
//
//  T time{0};
//  T increment{0};
//  T incrementMax{0};
//
//  void merge(const Timepoint<T> &t) {
//    assert(t.time == time);
//
//    increment += t.increment;
//    incrementMax += t.incrementMax;
//  }
//
//  std::ostream &display(std::ostream &os) const {
//    os << "(t=" << time << "|∂=" << increment << "|∂max=" << incrementMax
//       << ")";
//    return os;
//  }
//};


//template <typename T = int> struct Datapoint {
//
//  Datapoint() {}
//  Datapoint(const T overflow, const T consumption, const T overlap,
//            const T slackUnder, const T available)
//      : overflow(overflow), consumption(consumption), overlap(overlap),
//        slackUnder(slackUnder), available(available) {}
//
//  T overflow{0};
//  T consumption{0};
//  T overlap{0};
//  T slackUnder{0};
//  T available{0};
//
//  void reset() {
//    overflow = 0;
//    consumption = 0;
//    overlap = 0;
//    slackUnder = 0;
//    available = 0;
//  }
//
//  bool empty() {
//    return (overflow == 0 and consumption == 0 and overlap == 0 and
//            slackUnder == 0 and available == 0);
//  }
//
//  std::ostream &display(std::ostream &os) const {
//    if (overflow > 0)
//      os << "overflow=" << overflow;
//    if (consumption > 0)
//      os << " consumption=" << consumption;
//    if (overlap > 0)
//      os << " overlap=" << overlap;
//    if (slackUnder > 0)
//      os << " slackUnder=" << slackUnder;
//    if (available > 0)
//      os << " available=" << available;
//
//    return os;
//  }
//};
//
//template <typename T>
//std::ostream &operator<<(std::ostream &os, const Timepoint<T> &x) {
//  return x.display(os);
//}
//
//template <typename T>
//std::ostream &operator<<(std::ostream &os, const Datapoint<T> &x) {
//  return x.display(os);
//}

//template<typename T>
//class Solver;

template <typename T> class CumulativeOverlapFinding : public Constraint<T> {
private:
  Solver<T> &solver;
  NumericVar<T> capacity;
  Interval<T> schedule;
  std::vector<Interval<T>> task;
  std::vector<NumericVar<T>> demand;

  Matrix<Literal<T>> not_before;

  bool sign{bound::lower};
//  bool timetable_reasoning{true};

  // mapping from the tasks that need to be adjusted to the corresponding
  // left-cut interval
  std::vector<int> prec;

  // mapping the lower bound of the set explaining the adjustment
  std::vector<int> contact;

  // tasks that have been found to be in conflict with a left-cut interval
  //    SparseSet<> in_conflict;
  std::vector<int> in_conflict;

  List<Timepoint<T>> profile;
  std::vector<Datapoint<T>> data;
  std::vector<int>
      est_; // pointer to the profile element used for est of the task
  std::vector<int>
      lst_; // pointer to the profile element used for lst of the task
  std::vector<int>
      ect_; // pointer to the profile element used for ect of the task
  std::vector<int>
      lct_; // pointer to the profile element used for lct of the task

  std::vector<int> est_shared; // points to the possibly shared element whose
  // t is the est of the task
  std::vector<int> lst_shared; // points to the possibly shared element whose
  // t is the lst of the task
  std::vector<int> ect_shared; // points to the possibly shared element whose
  // t is the ect of the task
  std::vector<int> lct_shared; // points to the possibly shared element whose
  // t is the lct of the task

  int leftcut_pointer;

  int &get_shared_pointer(const int evt) {
    int n{3};
    auto evti{evt - 1};
    auto i{(evti) / n};
    switch ((evti % n)) {
    case 0:
      return est_shared[i];
    case 1:
      return ect_shared[i];
    default:
      return lct_shared[i];
//    case 2:
//      return lct_shared[i];
//    default:
//      return lst_shared[i];
    }
  }

  std::vector<int> event_ordering;
  std::vector<int> lct_order;
  int alpha;
  int beta;

  // tools for adjustments
  std::vector<int>
      minEct; // the minimum ect among tasks of the set use for explanation
  std::vector<int> maxOverflow; // the number of energy units to schedule on the
  // upper part of the profile
  std::vector<bool> isFeasible; // use to check if the correctonding scheduling
  // is feasible (without preemption)

  std::vector<std::vector<Literal<T>>> explanation;
  Reversible<size_t> num_explanations;
  std::vector<Literal<T>> pruning;

#ifdef STATS
    long unsigned int num_prop{0};
    long unsigned int num_useful{0};
    long unsigned int num_pruning{0};
#endif

public:
  template <typename ItTask, typename ItNVar>
  CumulativeOverlapFinding(Solver<T> &solver, const Interval<T> sched,
                           const NumericVar<T> capacity, const ItTask beg_task,
                           const ItTask end_task, const ItNVar beg_dem,
                           Matrix<Literal<T>> lits
                           //, const bool tt
    );
  virtual ~CumulativeOverlapFinding();

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;
  T mindemand(const unsigned i) const;
  T maxdemand(const unsigned i) const;
  T minenergy(const unsigned i) const;
  bool hasFixedPart(const unsigned i) const;
//  bool isLstWithoutFixedPart(const unsigned i) const;
  //    T overlapedFixedPartEnergy(const unsigned i, const unsigned j) const;
  //    std::vector<T> overlapedFixedPartEnergy();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void forwardpropagate();

//  void adjustment();
//  void detection();
  void addPrime(const int i, const int j);
  void rmPrime(const int i, const int j);

  int add_j_representation(const int i, const int j);
  void rm_j_representation(const int i, const int j, const int ej);
  int add_i_representation(const int i, const int j);
  void rm_i_representation(const int i, const int ei);

  void addExternalFixedPart(const int i, const int j);
  void rmExternalFixedPart(const int i, const int j);

  T scheduleOmega(const int i, const T max_lct, const bool adjustment = false);
  int makeNewEvent(const T t);

  void computeBound(const int i, std::vector<T> externalEnergy);
  void doPruning();

  // initialise the profile with every tasks
  void initialiseProfile();
  void clearData();

  // set the current profile to contain the tasks lct_order[0],..,lct_order[i]
  void setLeftCutToTime(const T t);
  void shrinkLeftCutToTime(const T t);
  void growLeftCutToTime(const T t);
  void rmTask(const int i);
  void addTask(const int i);

  std::vector<T> energyFixedPartExternalTasks();
  T fixedPartEnergy(const int i, const int j);

  // function used in explanation
  void computeExplanation(const int i, const int j);
  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;

#ifdef DBG_COF
  void verify(const char *msg);
  int debug_flag{5}; // 0 nothing but pruning // 1 summary // 2 task processing
                     // // 3 profile // 4 setleftcut
#endif
};

template <typename T>
std::string CumulativeOverlapFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << "t" << std::left << std::setw(3) << task[i].id() << ": [" << est(i)
     << ".." << lct(i) << "] (" << mindemand(i) << "x" << minduration(i) << ")";
  return ss.str();
}

template <typename T>
std::string CumulativeOverlapFinding<T>::asciiArt(const int i) const {
  std::stringstream ss;
  ss << std::setw(3) << std::right << mindemand(i) << "x" << std::setw(3)
     << std::left << minduration(i) << " " << std::right;
  for (auto k{0}; k < est(i); ++k) {
    ss << " ";
  }
    auto est_i{est(i)};
    auto ect_i{ect(i)};
    
    if (est_i == -Constant::Infinity<T>) {
        ss << "..." ;
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

template <typename T> T CumulativeOverlapFinding<T>::est(const unsigned i) const {
  return (sign == bound::lower ? task[i].getEarliestStart(solver)
                               : -task[i].getLatestEnd(solver));
}

template <typename T> T CumulativeOverlapFinding<T>::lst(const unsigned i) const {
  return (sign == bound::lower ? task[i].getLatestStart(solver)
                               : -task[i].getEarliestEnd(solver));
}

template <typename T> T CumulativeOverlapFinding<T>::ect(const unsigned i) const {
  return (sign == bound::lower ? task[i].getEarliestEnd(solver)
                               : -task[i].getLatestStart(solver));
}

template <typename T> T CumulativeOverlapFinding<T>::lct(const unsigned i) const {
  return (sign == bound::lower ? task[i].getLatestEnd(solver)
                               : -task[i].getEarliestStart(solver));
}

template <typename T>
T CumulativeOverlapFinding<T>::minduration(const unsigned i) const {
  return task[i].minDuration(solver);
}

template <typename T>
T CumulativeOverlapFinding<T>::maxduration(const unsigned i) const {
  return task[i].maxDuration(solver);
}

template <typename T>
T CumulativeOverlapFinding<T>::mindemand(const unsigned i) const {
  return demand[i].min(solver);
}

template <typename T>
T CumulativeOverlapFinding<T>::maxdemand(const unsigned i) const {
  return demand[i].max(solver);
}

template <typename T>
T CumulativeOverlapFinding<T>::minenergy(const unsigned i) const {
  return mindemand(i) * minduration(i);
}

template <typename T>
bool CumulativeOverlapFinding<T>::hasFixedPart(const unsigned i) const {
  return lst(i) < ect(i);
}

template <typename T>
template <typename ItTask, typename ItNVar>
CumulativeOverlapFinding<T>::CumulativeOverlapFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem,
    Matrix<Literal<T>> lits
                                                      //, const bool tt
)
    : solver(solver), not_before(std::move(lits))
//, timetable_reasoning(tt)
     , num_explanations(0, &(solver.getEnv())) {
  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    task.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  auto ip{task.size()};

  prec.resize(ip + 1, -1);
  est_.resize(ip);
  lst_.resize(ip);
  ect_.resize(ip);
  lct_.resize(ip);

  est_shared.resize(ip);
  lst_shared.resize(ip);
  ect_shared.resize(ip);
  lct_shared.resize(ip);

  contact.resize(ip + 1, -1);
  minEct.resize(ip, 0);
  isFeasible.resize(ip, false);
  maxOverflow.resize(ip, 0);

  for (unsigned i = 0; i < ip; ++i) {
    est_[i] = profile.create_element();
    ect_[i] = profile.create_element();
    lct_[i] = profile.create_element();
//    if (timetable_reasoning) {
//      lst_[i] = profile.create_element();
//    }
  }
  data.resize(profile.size() + 1);

  lct_order.resize(task.size());
  std::iota(lct_order.begin(), lct_order.end(), 0);
  event_ordering.resize(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);

  //        std::cout << capacity.max(solver) << std::endl;
  //        for(unsigned i{0}; i<task.size(); ++i) {
  //
  //            std::cout << std::endl << std::endl;
  //            for(unsigned j{0}; j<task.size(); ++j) {
  //                std::cout << task[i].start.id() << ":" <<
  //                demand[i].min(solver) << " <> " << task[j].start.id() << ":"
  //                << demand[j].min(solver)
  //                 << " -- " << solver.pretty(not_before(i,j)) << " // "
  //                 << solver.pretty(~not_before(i,j)) << std::endl;
  //            }
  //
  //        }
  //
  //        exit(1);
}

template <typename T> CumulativeOverlapFinding<T>::~CumulativeOverlapFinding() {}

template <typename T> void CumulativeOverlapFinding<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_CEDGEFINDING
  if (DBG_CEDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < task.size(); ++i) {
    solver.wake_me_on(lb<T>(task[i].end.id()), this->id());
    solver.wake_me_on(ub<T>(task[i].end.id()), this->id());
    solver.wake_me_on(lb<T>(demand[i].id()), this->id());
  }
}

template <typename T>
bool CumulativeOverlapFinding<T>::notify(const Literal<T>, const int) {
  return true;
}

template <typename T> void CumulativeOverlapFinding<T>::clearData() {
  auto p{profile.begin()};
  auto e{profile.end()};
  while (p != e) {
    data[p.index].reset();
    ++p;
  }
}

//template <typename T> bool CumulativeOverlapFinding<T>::isLstWithoutFixedPart(const unsigned evt) const {
//    if(timetable_reasoning and ((evt-1) % 4) == 3) {
//        if(not hasFixedPart((evt-1) / 4)) {
//            return true;
//        }
//    }
//    return false;
//}

template <typename T> void CumulativeOverlapFinding<T>::initialiseProfile() {

  profile.clear();

  // initialise the timepoints with the values from the domains
  for (auto i : lct_order) {
    est_shared[i] = est_[i];
    ect_shared[i] = ect_[i];
    lct_shared[i] = lct_[i];

    profile[est_[i]].time = est(i);
    profile[ect_[i]].time = ect(i);
    profile[lct_[i]].time = lct(i);

    profile[est_[i]].increment = profile[est_[i]].incrementMax = 0;
    profile[ect_[i]].increment = profile[ect_[i]].incrementMax = 0;
    profile[lct_[i]].increment = profile[lct_[i]].incrementMax = 0;

//    if (timetable_reasoning) {
//      lst_shared[i] = lst_[i];
//      profile[lst_[i]].time = lst(i);
//      profile[lst_[i]].increment = profile[lst_[i]].incrementMax = 0;
//    }
  }

  std::sort(event_ordering.begin(), event_ordering.end(),
            [this](const int i, const int j) {
              return this->profile[i].time < this->profile[j].time;
            });

  // populate the profile and merges duplicate times
  profile.add_front(event_ordering[0]);
  auto previous{event_ordering[0]};

  for (unsigned i{1}; i < event_ordering.size(); ++i) {
    auto current{event_ordering[i]};
//    if (isLstWithoutFixedPart(current))
//      continue;

    if (profile[previous].time == profile[current].time) {
      get_shared_pointer(current) = previous;
    } else {
      get_shared_pointer(current) = current;
      profile.add_after(previous, current);
      previous = current;
    }
  }

  // add tasks in lct order and make overload checks along the way
  leftcut_pointer = 0;

}

template <typename T> void CumulativeOverlapFinding<T>::rmTask(const int i) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "RM " << task[i].id() << std::endl;
  }
#endif

  profile[est_shared[i]].increment -= mindemand(i);
  profile[est_shared[i]].incrementMax -= mindemand(i);
  profile[ect_shared[i]].increment += mindemand(i);
  profile[lct_shared[i]].incrementMax += mindemand(i);
}

template <typename T> void CumulativeOverlapFinding<T>::addTask(const int i) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "ADD " << task[i].id() << std::endl;
  }
#endif

  profile[est_shared[i]].increment += mindemand(i);
  profile[est_shared[i]].incrementMax += mindemand(i);
  profile[ect_shared[i]].increment -= mindemand(i);
  profile[lct_shared[i]].incrementMax -= mindemand(i);
}

template <typename T>
void CumulativeOverlapFinding<T>::setLeftCutToTime(const T t) {
    
#ifdef DBG_COF
    if (DBG_COF and debug_flag > 4) {
        std::cout << "[set left cut to time " << t << "]\n";
    }
#endif
    
    if(t < -Constant::Infinity<T>) {
        int n{static_cast<int>(task.size())};
        while (leftcut_pointer < n) {
          addTask(lct_order[leftcut_pointer++]);
        }
        
#ifdef DBG_COF
        if (DBG_COF and debug_flag > 4) {
          std::cout << "[grown full (to " << leftcut_pointer << "-th task)]\n";
        }
        verify("after grow-full");
#endif

    } else if (leftcut_pointer < static_cast<int>(lct_order.size()) and
               lct(lct_order[leftcut_pointer]) < t) {
      growLeftCutToTime(t);
    } else {
      shrinkLeftCutToTime(t);
    }
}

template <typename T> void CumulativeOverlapFinding<T>::shrinkLeftCutToTime(const T t) {

#ifdef DBG_COF
  verify("before shrink-leftcut");
  if (DBG_COF and debug_flag > 4) {
    std::cout << "[shrink (cur=" << leftcut_pointer << ")]\n";
    }
#endif
    // remove the tasks that have a lct equal to or larger than t
    while (leftcut_pointer-- > 0 and t <= lct(lct_order[leftcut_pointer])) {
      rmTask(lct_order[leftcut_pointer]);
    }
    ++leftcut_pointer;

#ifdef DBG_COF
    if (DBG_COF and debug_flag > 4) {
        std::cout << "[shrank to " << leftcut_pointer << "-th task]\n";
    }
    verify("after shrink-leftcut");
#endif
}

template <typename T> void CumulativeOverlapFinding<T>::growLeftCutToTime(const T t) {

#ifdef DBG_COF
  verify("before grow-leftcut");
  if (DBG_COF and debug_flag > 4) {
    std::cout << "[grow (cur=" << leftcut_pointer << ")]\n";
    }
#endif

    // add the tasks that have a lct strictly smaller than t
    int n{static_cast<int>(task.size())};
    while (leftcut_pointer < n and t > lct(lct_order[leftcut_pointer])) {
      addTask(lct_order[leftcut_pointer++]);
    }

#ifdef DBG_COF
    if (DBG_COF and debug_flag > 4) {
        std::cout << "[grown to " << leftcut_pointer << "-th task]\n";
    }
    verify("after grow-leftcut");
#endif
}

template <typename T> void CumulativeOverlapFinding<T>::propagate() {

  sign = bound::lower;

  do {
    pruning.clear();

    std::sort(lct_order.begin(), lct_order.end(),
              [this](const int i, const int j) { return lct(i) < lct(j); });

#ifdef DBG_COF
    if (DBG_COF and debug_flag > 0) {
      std::cout << "\n\nstart ("
                << (sign == bound::lower ? "forward" : "backward")
                << ") COF propagation (" << capacity.max(solver) << ")\n";
      for (auto j : lct_order) {
        std::cout << "task " << std::setw(3) << task[j].id() << ": "
                  << asciiArt(j) << std::endl;
      }
    }
#endif

    initialiseProfile();

#ifdef DBG_COF
    verify("after init-profile");
#endif

    forwardpropagate();
    break;

    doPruning();
    sign ^= 1;
  } while (sign != bound::lower);
}

template <typename T> void CumulativeOverlapFinding<T>::forwardpropagate() {
    for (auto j : lct_order) {
        
#ifdef DBG_COF
      if (DBG_COF and debug_flag > 0) {
        std::cout << "* explore task " << std::left << std::setw(3)
                  << task[j].id() << std::endl;
      }
#endif
        
        setLeftCutToTime(lct(j) + Gap<T>::epsilon());
        auto ect_omega{scheduleOmega(j, lct(j))};
        
#ifdef DBG_COF
        if (DBG_COF and debug_flag > 0) {
          std::cout << "* ect(omega) = " << ect_omega << " / " << lct(j)
                    << std::endl;
        }
#endif
        
        if (ect_omega > lct(j)) {
            auto h{static_cast<hint>(num_explanations)};
            computeExplanation(-1, j);
            throw Failure<T>({this, h});
        } else {
            
            rmTask(j);
            
            for (auto i : lct_order) {
                if (i != j and solver.boolean.isUndefined(not_before(i,j).variable()) and (ect(i) > est(j) or lct(i) > lst(j))) {
                    
#ifdef DBG_COF
                  if (DBG_COF and debug_flag > 0) {
                    std::cout << " - assume task " << std::left << std::setw(3)
                              << task[i].id() << " << task " << task[j].id()
                              << std::endl;
                  }
#endif
                    
                    auto saved_size{profile.size()};
                    
                    auto ect_j_event{add_j_representation(i,j)};

                    if (lct(i) <= lct(j)) {
                        rmTask(i);
                    }
                    auto lct_i_event{add_i_representation(i,j)};

                    auto ect_h{scheduleOmega(j, lct(j))};
                    if (ect_h > lct(j)) {
                      // filter the variable b_i_j to false
                      //                        std::cout << "pruning!!: " <<
                      //                        solver.pretty(not_before(i,j))
                      //                        << "\n";
                      //                        solver.set(not_before(i,j),
                      //                        {this})
                      pruning.push_back(not_before(i, j));

                      auto h{num_explanations};
                      computeExplanation(i,j);

                      if (ect(i) > est(j)) {
                        explanation[h].push_back(task[i].start.after(est(i)));
                      }
                      if (lct(i) > lst(j)) {
                        explanation[h].push_back(task[j].end.before(lct(j)));
                      }

                      //                        exit(1);
                    }

                    //                    else {
                    //                        std::cout << "no pruning\n";
                    //                    }
                    rm_i_representation(i,lct_i_event);
                    rm_j_representation(i,j,ect_j_event);
                    
                    
                    if (lct(i) <= lct(j)) {
                        addTask(i);
                    }
                    
                    
                    while (profile.size() > saved_size) {
                      profile.pop_back();
                    }

                    //                    std::cout << "profile after:\n" <<
                    //                    profile << std::endl;
                }
            }
            
            addTask(j);
        }
    }
}

template <typename T> int CumulativeOverlapFinding<T>::add_j_representation(const int i, const int j) {
    
#ifdef DBG_COF
    if (DBG_COF and debug_flag > 0) {
      std::cout << "  - add rep for " << task[j].id() << std::endl;
    }
#endif
    
    int new_ect_j_event{-1};
    profile[ect_shared[i]].increment += mindemand(j);
    profile[ect_shared[i]].incrementMax += mindemand(j);
    auto t{profile.at(ect_shared[j])};
    auto target{ect(i) + minduration(j)};

    //    std::cout << " ect(j) should be " << target << std::endl;

    while (t->time < target)
        ++t;
    if (t->time == target) {

      //        std::cout << " found it! " << std::endl;

      t->increment -= mindemand(j);
      new_ect_j_event = t.index;
    } else {
        
        
        
        new_ect_j_event = makeNewEvent(target);
        --t;
        profile.add_after(t.index, new_ect_j_event);
        profile[new_ect_j_event].increment -= mindemand(j);

        //        std::cout << " must create it " << profile[new_ect_j_event] <<
        //        " and insert after " << profile[t.index] << std::endl;
    }
    profile[lct_shared[j]].incrementMax -= mindemand(j);
    return new_ect_j_event;
}

template <typename T> void CumulativeOverlapFinding<T>::rm_j_representation(const int i, const int j, const int new_ect_j_event) {
    profile[ect_shared[i]].increment -= mindemand(j);
    profile[ect_shared[i]].incrementMax -= mindemand(j);
    profile[new_ect_j_event].increment += mindemand(j);
    profile[lct_shared[j]].incrementMax += mindemand(j);
}

template <typename T> int CumulativeOverlapFinding<T>::add_i_representation(const int i, const int j) {
    int new_lct_i_event{-1};
    profile[est_shared[i]].increment += mindemand(i);
    profile[est_shared[i]].incrementMax += mindemand(i);
    profile[ect_shared[i]].increment -= mindemand(i);
    if(lct(i) <= lst(j)) {
        profile[lct_shared[i]].incrementMax -= mindemand(i);
        new_lct_i_event = lct_shared[i];
    } /*else if(timetable_reasoning) {
       profile[lst_shared[j]].incrementMax -= mindemand(i);
       new_lct_i_event = lst_shared[j];
       }*/
    else {
      auto target{lst(j)};
      auto t{profile.at(lct_shared[i])};
      while (target < t->time) {
        --t;
      }
      if (t->time == target) {
        t->incrementMax -= mindemand(i);
        new_lct_i_event = t.index;
      } else {
        new_lct_i_event = makeNewEvent(lst(j));
        profile.add_after(t.index, new_lct_i_event);
        profile[new_lct_i_event].incrementMax -= mindemand(i);
      }
    }
    return new_lct_i_event;
}

template <typename T> void CumulativeOverlapFinding<T>::rm_i_representation(const int i, const int new_lct_i_event) {
    profile[est_shared[i]].increment -= mindemand(i);
    profile[est_shared[i]].incrementMax -= mindemand(i);
    profile[ect_shared[i]].increment += mindemand(i);
    profile[new_lct_i_event].incrementMax += mindemand(i);
}


template <typename T> void CumulativeOverlapFinding<T>::addPrime(const int i, const int j) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "ADD " << task[i].id() << "'" << std::endl;
  }
#endif

  profile[est_shared[i]].increment += mindemand(i);
  profile[est_shared[i]].incrementMax += mindemand(i);

  if (ect(i) < lct(j)) {
    profile[ect_shared[i]].increment -= mindemand(i);
    profile[ect_shared[i]].incrementMax -= mindemand(i);
  } else {
    profile[lct_shared[j]].increment -= mindemand(i);
    profile[lct_shared[j]].incrementMax -= mindemand(i);
  }
}

template <typename T> void CumulativeOverlapFinding<T>::rmPrime(const int i, const int j) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "RM " << task[i].id() << "'" << std::endl;
  }
#endif

  profile[est_shared[i]].increment -= mindemand(i);
  profile[est_shared[i]].incrementMax -= mindemand(i);

  if (ect(i) < lct(j)) {
    profile[ect_shared[i]].increment += mindemand(i);
    profile[ect_shared[i]].incrementMax += mindemand(i);
  } else {
    profile[lct_shared[j]].increment += mindemand(i);
    profile[lct_shared[j]].incrementMax += mindemand(i);
  }
}


template <typename T> void CumulativeOverlapFinding<T>::addExternalFixedPart(const int i, const int j) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "ADD " << task[i].id() << "'" << std::endl;
  }
#endif
  for (auto k : lct_order) {
    if (lct(k) > lct(j) and lst(k) < ect(k) and lst(k) < lct(j) and k != i) {
      profile[lst_shared[k]].increment += mindemand(k);
      profile[lst_shared[k]].incrementMax += mindemand(k);

      if (ect(k) < lct(j)) {
        profile[ect_shared[k]].increment -= mindemand(k);
        profile[ect_shared[k]].incrementMax -= mindemand(k);
      } else {
        profile[lct_shared[j]].increment -= mindemand(k);
        profile[lct_shared[j]].incrementMax -= mindemand(k);
      }
    }
  }
}



template <typename T> void CumulativeOverlapFinding<T>::rmExternalFixedPart(const int i, const int j) {

#ifdef DBG_COF
  if (DBG_COF and debug_flag > 3) {
    std::cout << "ADD " << task[i].id() << "'" << std::endl;
  }
#endif
  for (auto k : lct_order) {
    if (lct(k) > lct(j) and lst(k) < ect(k) and lst(k) < lct(j) and k != i) {
      profile[lst_shared[k]].increment -= mindemand(k);
      profile[lst_shared[k]].incrementMax -= mindemand(k);

      if (ect(k) < lct(j)) {
        profile[ect_shared[k]].increment += mindemand(k);
        profile[ect_shared[k]].incrementMax += mindemand(k);
      } else {
        profile[lct_shared[j]].increment += mindemand(k);
        profile[lct_shared[j]].incrementMax += mindemand(k);
      }
    }
  }
}

template <typename T> T CumulativeOverlapFinding<T>::fixedPartEnergy(const int i, const int j){
    T E = 0;
    if (lct(i) > lct(j) and lst(i) < ect(i) and lst(i) < lct(j))
        E = std::max(0, std::min(ect(i), lct(j)) - lst(i)) * mindemand(i);
    return E;
}

template <typename T> std::vector<T> CumulativeOverlapFinding<T>::energyFixedPartExternalTasks(){
    std::vector<T> energy;
    int n = static_cast<int>(lct_order.size());
    energy.resize(n, 0);
    for (auto k{0}; k < n; ++k){
        auto j{lct_order[k]};
        int E = 0;
        for (auto l{k+1}; l < n; ++l){
            auto i{lct_order[l]};
            E += fixedPartEnergy(i,j);
        }
        energy[j] = E;
        while (k+1 < n and lct(lct_order[k+1]) == lct(lct_order[k])){
            auto q{lct_order[k+1]};
            energy[q] = energy[lct_order[k]];
            ++k;
        }
    }
    return energy;
}



template <typename T>
void CumulativeOverlapFinding<T>::computeBound(const int i, std::vector<T> externalEnergy) {
  T E{0};
  alpha = -1;
  beta = -1;
  T minSlack[2] = {Constant::Infinity<T>, Constant::Infinity<T>};
  // if (timetable_reasoning){
  for (auto j : lct_order) {
    if (lct(j) == lct(i))
      break;
    E += minenergy(j);
    if (lct(j) <= ect(i) and est(i) < lct(j)) {
      auto slack{(capacity.max(solver) - mindemand(i)) * lct(j) - E -
                 externalEnergy[j] + fixedPartEnergy(i, j)};
      if (slack < minSlack[0] and data[lct_shared[j]].overflow > 0) {
        minSlack[0] = slack;
        alpha = j;
      }
    } else if (lct(j) > ect(i)) {
      auto slack{capacity.max(solver) * lct(j) - E - externalEnergy[j] +
                 fixedPartEnergy(i, j)};
      if (slack < minSlack[1] and data[lct_shared[j]].overflow > 0) {
        minSlack[1] = slack;
        beta = j;
      }
    }
    }
}

template <typename T>
int CumulativeOverlapFinding<T>::makeNewEvent(const T t) {
  auto new_event{profile.create_element(t, 0, 0)};
  if (data.size() <= static_cast<size_t>(new_event)) {
    data.resize(new_event + 1);
  } else {
    data[new_event].reset();
  }
  return new_event;
}

template <typename T>
T CumulativeOverlapFinding<T>::scheduleOmega(const int i, const T max_lct,
                                             const bool adjustment) {

  auto saved_size{profile.size()};

#ifdef DBG_COF
    if (DBG_COF and debug_flag > 2) {
      std::cout << "[schedule tasks until t=" << max_lct
                << "] profile=" << std::endl
                << profile;
    }
#endif

    if (adjustment) {
      clearData();

#ifdef DBG_COF
      if (DBG_COF and debug_flag > 2) {
        auto next{profile.begin()};
        auto stop{profile.end()};

        while (next != stop and next->time <= max_lct) {
          auto t{next};
          ++next;
          std::cout << *t << ": " << data[t.index] << std::endl;
        }
      }
#endif
    }

    auto next{profile.begin()};

    T C{capacity.max(solver)};
    T overflow = 0;
    T omega_ect{-Constant::Infinity<T>};
    T S{0};
    T h_req{0};

    T overlap{0};
    T slackUnder{0};
    T available{0};
    T h{mindemand(i)};
    isFeasible[i] = true;

    auto stop{profile.end()};

    while (next != stop and next->time <= max_lct) {
      auto t{next};
      ++next;

#ifdef DBG_COF
      if (DBG_COF and debug_flag > 2) {
        std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
      }
#endif

      data[t.index].overflow = overflow;

      auto l = next->time - t->time;

      if (adjustment) {
        data[t.index].overlap = overlap;
        data[t.index].slackUnder = slackUnder;
        data[t.index].available = available;
      }

      // S is the sum of demands of the tasks that could be processed at time
      S += t->incrementMax;
      // h_max is the min between the resource's capacity and the total of the
      // demands
      auto h_max{std::min(S, C)};
      // h_req is the total demand counting tasks processed at their earliest
      h_req += t->increment;
      // h_cons is the amount of resource actually used in the optimistic
      // scenario (min between what is required + due from earlier, and what is
      // available)
      auto h_cons{std::min(h_req + overflow, h_max)};

      if (overflow > 0 and (h_cons - h_req) > 0) {
        auto al{std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req))};
        // there is some overflow, and it will be resorbed by the next time
        // point
        if (al < l) {
          l = al;
          auto new_event{makeNewEvent(t->time + l)};
          profile.add_after(t.index, new_event);

#ifdef DBG_COF
          if (DBG_COF and debug_flag > 4) {
            std::cout << " [create tmp event " << *(profile.at(new_event))
                      << " before " << *next << " b/c " << al << "] ";
          }
#endif

          next = profile.at(new_event);
        }
      }

      // overflow is the deficit on resource for that period (because tasks are
      // set to their earliest)profile[est_[ip]].time = _est;
      overflow += (h_req - h_cons) * l;

      if (overflow > 0)
        isFeasible[i] = false;

      if (adjustment) {
        overlap += (std::max(h_cons - (C - h), 0) * l);
        slackUnder += (std::max(std::min(C - h, h_max) - h_cons, 0) * l);
        available += (std::min(C - h_cons, h) * l);

        data[t.index].consumption = h_cons;
      }

      if (overflow > 0)
        omega_ect = Constant::Infinity<T>;
      else if (h_cons > 0)
        omega_ect = profile[profile.next(t.index)].time;

#ifdef DBG_COF
      if (DBG_COF and debug_flag > 2) {
        std::cout << " h_max=" << h_max << ", h_req=" << h_req
                  << ", h_cons=" << h_cons << " (" << data[t.index] << ")";
        if (omega_ect == Constant::Infinity<T>)
          std::cout << ", ect=inf";
        else if (omega_ect != -Constant::Infinity<T>)
          std::cout << ", ect=" << omega_ect;

        std::cout << std::endl;
      }

      assert(h_cons >= 0 or next->time <= max_lct);
#endif
    }

    if (not adjustment) {
      contact[i] = -1;

#ifdef DBG_COF
      if (DBG_COF and debug_flag > 4) {
        std::cout << "search contact\n";
      }
#endif

      auto t{next};
      if (overflow == 0 and omega_ect > schedule.end.min(solver)) {

#ifdef DBG_COF
        if (DBG_COF and debug_flag > 5) {
          std::cout << "go to overflow (bound pruning)\n";
        }
#endif

        auto limit{profile.begin()};
        while (data[t.index].overflow == 0 and t != limit) {
          --t;
        }
        overflow = data[t.index].overflow;
#ifdef DBG_COF
        if (DBG_COF and debug_flag > 5) {
          std::cout << "hack overflow " << *t << ": " << data[t.index]
                    << std::endl;
        }
#endif
      }

      if (overflow > 0) {
        contact[i] = profile.begin().index;

#ifdef DBG_COF
        if (DBG_COF and debug_flag > 5) {
          std::cout << "normal contact\n";
        }
#endif

        do {
          --t;
          if (data[t.index].overflow < overflow) {
            contact[i] = t.index;
            break;
          }
        } while (t != profile.begin());

#ifdef DBG_COF
        if (DBG_COF and debug_flag > 2) {
          std::cout << "contact[" << task[i].id()
                    << "] = " << profile[contact[i]] << " (" << contact[i]
                    << ") contact time = " << profile[contact[i]].time << ": "
                    << data[contact[i]] << std::endl;
        }
#endif
      }
    }

    //    contact[i] = resetProfile(contact[i], saved_size);
    while (profile.size() > saved_size) {
      profile.pop_back();
    }

    return omega_ect;
}

template <typename T> void CumulativeOverlapFinding<T>::doPruning() {

#ifdef STATS
    num_pruning += pruning.size();
    if(not pruning.empty()) {
        ++num_useful;
        //        std::cout << "EF" << this->id() << ": " << num_useful << "/"
        //        << num_prop << " (" << num_pruning << ")\n"; std::cout << "EF"
        //        << this->id() << ": " <<
        //        static_cast<double>(num_useful)/static_cast<double>(num_prop)
        //        << " (" << num_pruning << ")\n";
    }
#endif

    auto h{static_cast<hint>(num_explanations - pruning.size())};
    for (auto p : pruning) {
      solver.set(p, {this, h++});
    }
}

template <typename T>
void CumulativeOverlapFinding<T>::computeExplanation(const int x, const int y) {
//  auto n{static_cast<int>(task.size())};

  auto h{num_explanations};
  if (explanation.size() <= h) {
    explanation.resize(h + 1);
  } else {
    explanation[h].clear();
  }
  ++num_explanations;
  auto k{leftcut_pointer - 1};
  auto j{lct_order[k]};


//  assert(i == n or lct(j) == lct(prec[i]));

  auto t{profile[contact[y]].time};

#ifdef DBG_COF
  if (DBG_COF and debug_flag >= 0) {
    if (x != -1) {
        std::cout << "*** adjustment task " << task[x].id() << " not before task " << task[y].id() << std::endl;
    } else {
      std::cout << "*** overload fail";
    }
    std::cout << " capacity = " << capacity.max(solver) << " @propag #"
              << solver.num_cons_propagations << "\n";
  }
#endif

    
    
  do {
    if (ect(j) > t and j != x and j != y) {

        auto saved_sign{sign};
        sign = bound::lower;
      explanation[h].push_back(task[j].start.after(est(j)));
      explanation[h].push_back(task[j].end.before(lct(j)));
      sign = saved_sign;

#ifdef DBG_COF
      if (DBG_COF and debug_flag >= 0) {
        std::cout << "task " << std::setw(3) << task[j].id() << ": "
                  << asciiArt(j) << std::endl;
      }
#endif
    }
  } while (k > 0 and lct(j = lct_order[--k]) > t);
    
    
    if(x == -1) {
        auto saved_sign{sign};
        sign = bound::lower;
     
        explanation[h].push_back(task[y].start.after(est(y)));
        explanation[h].push_back(task[y].end.before(lct(y)));
        
        sign = saved_sign;
    }

#ifdef DBG_COF
    if (DBG_COF and debug_flag >= 0) {
      std::cout << std::endl;
    }
#endif
}

template <typename T>
void CumulativeOverlapFinding<T>::xplain(const Literal<T>, const hint h,
                                         std::vector<Literal<T>> &Cl) {

  for (auto p : explanation[h]) {
    Cl.push_back(p);
  }
}

template <typename T>
std::ostream &CumulativeOverlapFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding";

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

template <typename T>
std::ostream &CumulativeOverlapFinding<T>::print_reason(std::ostream &os, const hint) const {
  os << "cumulative-edge-finding";
  return os;
}


#ifdef DBG_COF // debugging stuff
template <typename T>
void CumulativeOverlapFinding<T>::verify(const char* msg) {
  std::vector<Timepoint<T>> events;
  for (int k{0}; k < leftcut_pointer; ++k) {
    auto i{lct_order[k]};
    events.emplace_back(est(i), mindemand(i), mindemand(i));
    events.emplace_back(ect(i), -mindemand(i), 0);
    events.emplace_back(lct(i), 0, -mindemand(i));
  }
  std::sort(events.begin(), events.end(),
            [](const Timepoint<T> &a, const Timepoint<T> &b) {
              return a.time < b.time;
            });

  auto tp{profile.begin()};
  auto i{tp};
  while (i != profile.end()) {
    ++i;
    if (i == tp) {
      std::cout << "infinite loop " << msg << "\n"; //
      exit(1);
    }
  }

  auto tv{events.begin()};

  while (tv != events.end()) {
    T d{0};
    T dv{0};
    T dm{0};
    T dmv{0};

    auto now{tv->time};
    do {
      dv += tv->increment;
      dmv += tv->incrementMax;
    } while (++tv != events.end() and tv->time == now);

    bool bug{false};
    while (tp->time < now) {
      if (tp->increment != 0 or tp->incrementMax != 0) {
        std::cout << "non-zero non-task event!\n";
        bug = true;
      }
      ++tp;
    }

    d += tp->increment;
    dm += tp->incrementMax;

    if (d != dv) {
      std::cout << "discrepancy increment @" << *tp << "|" << *(tv - 1) << ": "
                << d << "/" << dv << std::endl;
      bug = true;
    }
    if (dm != dmv) {
      std::cout << "discrepancy increment max @" << *tp << "|" << *(tv - 1)
                << ": " << dm << "/" << dmv << std::endl;
      bug = true;
    }

    ++tp;

    if (bug) {
      std::cout << "bug " << msg << "\n"; //
      auto itp{profile.begin()};
      auto itv{events.begin()};

      while (true) {
        if (itp != tp) {
          while (itp->increment == 0 and itp->incrementMax == 0) {
            ++itp;
          }
        }
        if (itp != tp) {
          while (itv != tv and itv->time <= itp->time) {
            std::cout << *itp << " / " << *itv << std::endl;
            ++itv;
          }
        }
        ++itp;
        std::cout << std::endl;
        if (itp == tp and itv == tv)
          break;

        if (itp == tp or itv == tv) {
          std::cout << "rebug\n";
          exit(1);
        }
      }
      exit(1);
    }
  }
}
#endif


} // namespace tempo

#endif //CUMULATIVEOVERLAPFINDING_HPP
