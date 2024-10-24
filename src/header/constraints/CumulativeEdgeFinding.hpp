/************************************************
 * Tempo CumulativeEdgeFinding.hpp
 * Implementation of the "strong edge-finding" algorithm as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder
 *rule for the cumulative constraint. In Proceedings of the Twenty-Fifth
 *International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI
 *Press, 3103–3109.
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

#ifndef TEMPO_CUMULATIVEEDGEFINDING_HPP
#define TEMPO_CUMULATIVEEDGEFINDING_HPP

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

template <typename T = int> struct Timepoint {

  Timepoint() {}
  Timepoint(T time, T increment, T incrementMax)
      : time(time), increment(increment), incrementMax(incrementMax) {}

  T time{0};
  T increment{0};
  T incrementMax{0};

  void merge(const Timepoint<T> &t) {
    assert(t.time == time);

    increment += t.increment;
    incrementMax += t.incrementMax;
  }

  std::ostream &display(std::ostream &os) const {
    os << "(t=" << time << "|∂=" << increment << "|∂max=" << incrementMax
       << ")";
    return os;
  }
};


template <typename T = int> struct Datapoint {

  Datapoint() {}
  Datapoint(const T overflow, const T consumption, const T overlap,
            const T slackUnder, const T available)
      : overflow(overflow), consumption(consumption), overlap(overlap),
        slackUnder(slackUnder), available(available) {}

  T overflow{0};
  T consumption{0};
  T overlap{0};
  T slackUnder{0};
  T available{0};

  std::ostream &display(std::ostream &os) const {
    if (overflow > 0)
      os << "overflow=" << overflow;
    if (consumption > 0)
      os << " consumption=" << overflow;
    if (overlap > 0)
      os << " overlap=" << overlap;
    if (slackUnder > 0)
      os << " slackUnder=" << slackUnder;
    if (available > 0)
      os << " available=" << available;

    return os;
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Timepoint<T> &x) {
  return x.display(os);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Datapoint<T> &x) {
  return x.display(os);
}

template<typename T>
class Solver;

template <typename T> class CumulativeEdgeFinding : public Constraint<T> {
private:
  Solver<T> &solver;
  NumericVar<T> capacity;
  Interval<T> schedule;
  std::vector<Interval<T>> task;
  std::vector<NumericVar<T>> demand;

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
      ect_; // pointer to the profile element used for ect of the task
  std::vector<int>
      lct_; // pointer to the profile element used for lct of the task

  std::vector<int> est_shared; // points to the possibly shared element whose
  // t is the est of the task
  std::vector<int> ect_shared; // points to the possibly shared element whose
  // t is the ect of the task
  std::vector<int> lct_shared; // points to the possibly shared element whose
  // t is the lct of the task

  int leftcut_pointer;

  int &get_shared_pointer(const int evt) {
    auto evti{evt - 1};
    auto i{(evti) / 3};
    switch ((evti % 3)) {
    case 0:
      return est_shared[i];
    case 1:
      return ect_shared[i];
    default:
      return lct_shared[i];
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

public:
  template <typename ItTask, typename ItNVar>
  CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched,
                        const NumericVar<T> capacity, const ItTask beg_task,
                        const ItTask end_task, const ItNVar beg_dem);
  virtual ~CumulativeEdgeFinding();

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
  //    bool hasFixedPart(const unsigned i) const;
  //    T overlapedFixedPartEnergy(const unsigned i, const unsigned j) const;
  //    std::vector<T> overlapedFixedPartEnergy();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void forwardAdjustment();
  void forwardDetection();
  void addPrime(const int i, const int j);
  void rmPrime(const int i, const int j);
  T scheduleOmega(const int i, const T max_lct, const bool adjustment = false);
  int makeNewEvent(const T t);

  void computeBound(const int i);
  void doPruning();

  // initialise the profile with every tasks
  void initialiseProfile();

  // set the current profile to contain the tasks lct_order[0],..,lct_order[i]
  void setLeftCutToTime(const T t);
  void shrinkLeftCutToTime(const T t);
  void growLeftCutToTime(const T t);
  void rmTask(const int i);
  void addTask(const int i);

  // function used in explanation
   void computeForwardExplanation(const int i);
  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;
};

template <typename T>
std::string CumulativeEdgeFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << "t" << task[i].id() << ": [" << est(i) << ".." << lct(i) << "] ("
     << mindemand(i) << "x" << minduration(i) << ")";
  return ss.str();
}

template <typename T>
std::string CumulativeEdgeFinding<T>::asciiArt(const int i) const {
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

template <typename T> T CumulativeEdgeFinding<T>::est(const unsigned i) const {
  return task[i].getEarliestStart(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lst(const unsigned i) const {
  return task[i].getLatestStart(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::ect(const unsigned i) const {
  return task[i].getEarliestEnd(solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lct(const unsigned i) const {
  return task[i].getLatestEnd(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::minduration(const unsigned i) const {
  return task[i].minDuration(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxduration(const unsigned i) const {
  return task[i].maxDuration(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::mindemand(const unsigned i) const {
  return demand[i].min(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxdemand(const unsigned i) const {
  return demand[i].max(solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::minenergy(const unsigned i) const {
  return mindemand(i) * minduration(i);
}

//template <typename T>
//bool CumulativeEdgeFinding<T>::hasFixedPart(const unsigned i) const {
//  return lst(i) < ect(i);
//}

// template <typename T>
// T CumulativeEdgeFinding<T>::overlapedFixedPartEnergy(const unsigned i, const
// unsigned j) const {
//     auto energy{0};
//     if (lct(i) > lct(j) and hasFixedPart(i) and lst(i) < lct(j)) {
//         energy = (std::max(0, std::min(ect(i), lct(j)) - lst(i))) *
//         mindemand(i);
//     }
//     return energy;
// }

// template <typename T>
// std::vector<T> CumulativeEdgeFinding<T>::overlapedFixedPartEnergy() {
//     std::vector<T> energy;
//     energy.resize(task.size());
//     for (auto j : lct_order) {
//         auto E{0};
//         for (auto i = j + 1; i < task.size(); ++i) {
//             E += overlapedFixedPartEnergy(i, j);
//         }
//         energy[j] = E;
//         while (j + 1 < task.size() and lct(j + 1) == lct(j)) {
//             energy[j + 1] = energy[j];
//             ++j;
//         }
//     }
//     return energy;
// }

template <typename T>
template <typename ItTask, typename ItNVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem)
    : solver(solver), num_explanations(0, &(solver.getEnv())) {
  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    task.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  auto ip{task.size()};

  //    in_conflict.reserve(ip);
  prec.resize(ip);
  est_.resize(ip);
  ect_.resize(ip);
  lct_.resize(ip);

  est_shared.resize(ip);
  ect_shared.resize(ip);
  lct_shared.resize(ip);

  contact.resize(ip, -1);
  minEct.resize(ip, 0);
  isFeasible.resize(ip, false);
  maxOverflow.resize(ip, 0);

  for (unsigned i = 0; i < ip; ++i) {
    est_[i] = profile.create_element();
    ect_[i] = profile.create_element();
    lct_[i] = profile.create_element();
  }
  data.resize(profile.size() + 1);

  lct_order.resize(task.size());
  std::iota(lct_order.begin(), lct_order.end(), 0);
  event_ordering.resize(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}


template <typename T> void CumulativeEdgeFinding<T>::post(const int idx) {

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
bool CumulativeEdgeFinding<T>::notify(const Literal<T>, const int) {
  return true;
}

template <typename T> void CumulativeEdgeFinding<T>::initialiseProfile() {

  profile.clear();

  // initialise the timepoints with the values from the domains
  for (auto i : lct_order) {
    est_shared[i] = est_[i];
    ect_shared[i] = ect_[i];
    lct_shared[i] = lct_[i];

    profile[est_[i]].time = est(i);
    profile[ect_[i]].time = ect(i);
    profile[lct_[i]].time = lct(i);
    profile[est_[i]].increment = mindemand(i);
    profile[est_[i]].incrementMax = mindemand(i);
    profile[ect_[i]].increment = -mindemand(i);
    profile[lct_[i]].incrementMax = -mindemand(i);
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
    if (profile[previous].time == profile[current].time) {
      profile[previous].merge(profile[current]);
      get_shared_pointer(current) = previous;
    } else {
      get_shared_pointer(current) = current;
      profile.add_after(previous, current);
      previous = current;
    }
  }

//  for (auto i : lct_order) {
//    assert(profile[est_shared[i]].time == est(i));
//    assert(profile[ect_shared[i]].time == ect(i));
//    assert(profile[lct_shared[i]].time == lct(i));
//  }

  leftcut_pointer = static_cast<int>(lct_order.size());
}


template <typename T> void CumulativeEdgeFinding<T>::rmTask(const int i) {
  profile[est_shared[i]].increment -= mindemand(i);
  profile[est_shared[i]].incrementMax -= mindemand(i);
  profile[ect_shared[i]].increment += mindemand(i);
  profile[lct_shared[i]].incrementMax += mindemand(i);
}

template <typename T> void CumulativeEdgeFinding<T>::addTask(const int i) {
  profile[est_shared[i]].increment += mindemand(i);
  profile[est_shared[i]].incrementMax += mindemand(i);
  profile[ect_shared[i]].increment -= mindemand(i);
  profile[lct_shared[i]].incrementMax -= mindemand(i);
}

template <typename T> void CumulativeEdgeFinding<T>::setLeftCutToTime(const T t) {
  if (leftcut_pointer < static_cast<int>(lct_order.size()) and
      lct(lct_order[leftcut_pointer]) < t) {
    growLeftCutToTime(t);
  } else {
    shrinkLeftCutToTime(t);
  }
}

template <typename T> void CumulativeEdgeFinding<T>::shrinkLeftCutToTime(const T t) {

//    std::cout << "shrink until " << t << ":\n" << profile << std::endl;

// remove the tasks that have a lct equal to or larger than t
while (leftcut_pointer-- > 0 and t <= lct(lct_order[leftcut_pointer])) {

  //        std::cout << " - remove " << prettyTask(lct_order[leftcut_pointer])
  //        << std::endl;

  rmTask(lct_order[leftcut_pointer]);
}

//    std::cout << "==>\n" << profile << std::endl;

++leftcut_pointer;
}

template <typename T> void CumulativeEdgeFinding<T>::growLeftCutToTime(const T t) {
  // add the tasks that have a lct strictly smaller than t
  while (t > lct(lct_order[leftcut_pointer])) {
    addTask(lct_order[leftcut_pointer++]);
  }
}

template <typename T> void CumulativeEdgeFinding<T>::propagate() {

std::sort(lct_order.begin(), lct_order.end(),
          [this](const int i, const int j) { return lct(i) < lct(j); });

#ifdef DBG_SEF
if (DBG_SEF) {
  std::cout << "\n\nstart propagation (" << capacity.max(solver) << ")\n";
  for (auto j : lct_order) {
    std::cout << "task " << task[j].id() << ": " << asciiArt(j) << std::endl;
  }
}
#endif

initialiseProfile();

forwardDetection();

forwardAdjustment();

doPruning();
}

template <typename T> void CumulativeEdgeFinding<T>::addPrime(const int i, const int j) {
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

template <typename T> void CumulativeEdgeFinding<T>::rmPrime(const int i, const int j) {
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

template <typename T>
void CumulativeEdgeFinding<T>::computeBound(const int i) {
  T E{0};
  alpha = -1;
  beta = -1;
  T minSlack[2] = {Constant::Infinity<T>, Constant::Infinity<T>};
  for (auto j : lct_order) {
    if (lct(j) == lct(i))
      break;
    E += minenergy(j);
    if (lct(j) <= ect(i) and est(i) < lct(j)) {
      auto slack{(capacity.max(solver) - mindemand(i)) * lct(j) - E};
      if (slack < minSlack[0] and data[lct_[j]].overflow > 0) {
        minSlack[0] = slack;
        alpha = j;
      }
    } else if (lct(j) > ect(i)) {
      auto slack{capacity.max(solver) * lct(j) - E};
      if (slack < minSlack[1] and data[lct_[j]].overflow > 0) {
        minSlack[1] = slack;
        beta = j;
      }
    }
  }
}

template <typename T>
int CumulativeEdgeFinding<T>::makeNewEvent(const T t) {
  auto new_event{profile.create_element(t, 0, 0)};
  if (data.size() <= static_cast<size_t>(new_event)) {
    data.resize(new_event + 1);
  }
  return new_event;
}

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmega(const int i, const T max_lct,
                                          const bool adjustment) {

  auto saved_size{profile.size()};

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "[schedule tasks until t=" << max_lct
              << "] profile=" << std::endl
              << profile;
  }
#endif

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

  while (next->time <= max_lct) {
    //     while (next != stop) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
    }
#endif
      
    data[t.index].overflow = overflow;

    auto l = next->time - t->time;

    // S is the sum of demands of the tasks that could be processed at time
    S += t->incrementMax;
    // h_max is the min between the resource's capacity and the total of the
    // demands
    auto h_max{std::min(S, C)};
    // h_req is the total demand counting tasks processed at their earliest
    h_req += t->increment;
    // h_cons is the amount of resource actually used in the optimistic scenario
    // (min between what is required + due from earlier, and what is available)
    auto h_cons{std::min(h_req + overflow, h_max)};

    if (overflow > 0 and (h_cons - h_req) != 0) {
      auto al{std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req))};
      // there is some overflow, and it will be resorbed by the next time point
      if (al < l) {
        l = al;
        auto new_event{makeNewEvent(t->time + l)};
        profile.add_after(t.index, new_event);
        next = profile.at(new_event);
      }
    }

    // overflow is the deficit on resource for that period (because tasks are
    // set to their earliest)profile[est_[ip]].time = _est;
    overflow += (h_req - h_cons) * l;

    if (overflow > 0)
      isFeasible[i] = false;

    
    if (adjustment) {
      data[t.index].consumption = h_cons;
      overlap += (std::max(h_cons - (C - h), 0) * l);
      slackUnder += (std::max(std::min(C - h, h_max) - h_cons, 0) * l);
      available += (std::min(C - h_cons, h) * l);

      data[t.index].overlap = overlap;
      data[t.index].slackUnder = slackUnder;
      data[t.index].available = available;
    }

    if (overflow > 0)
      omega_ect = Constant::Infinity<T>;
    else if (h_cons > 0)
      omega_ect = profile[profile.next(t.index)].time;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << " h_max=" << h_max << ", h_req=" << h_req
                << ", h_cons=" << h_cons << " (" << data[t.index] << ")";
      if (omega_ect == Constant::Infinity<T>)
        std::cout << ", ect=inf";
      else if (omega_ect != -Constant::Infinity<T>)
        std::cout << ", ect=" << omega_ect;

      std::cout << std::endl;
    }
#endif
  }

  if (not adjustment) {
    contact[i] = -1;
    if (overflow > 0) {
      contact[i] = profile.begin().index;
      auto t{next};
      do {
        --t;
        if (data[t.index].overflow < overflow) {
          contact[i] = t.index;
          break;
        }
      } while (t != profile.begin());

#ifdef DBG_SEF
      if (DBG_SEF) {
        std::cout << "contact[" << task[i].id() << "] = " << profile[contact[i]]
                  << " (" << contact[i]
                  << ") contact time = " << profile[contact[i]].time
                  << std::endl;
      }
#endif
    }
  }

  while (profile.size() > saved_size) {
    profile.pop_back();
  }

  return omega_ect;
}

template <typename T> void CumulativeEdgeFinding<T>::doPruning() {

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

template <typename T> void CumulativeEdgeFinding<T>::forwardAdjustment() {

#ifdef DBG_SEF
  if (DBG_SEF) {
    if (in_conflict.empty())
      std::cout << "\nno forward adjustment\n";
    else
      std::cout << "\nstart forward adjustment\n";
  }
#endif

  while (not in_conflict.empty()) {
    auto i{in_conflict.back()};
    in_conflict.pop_back();

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << " adjustment on " << task[i].id() << "\n";
    }
#endif

    auto j{prec[i]};

    // compute the profile without i, but to record the overlap and
    // slackUnder
    setLeftCutToTime(lct(j) + Gap<T>::epsilon());
    scheduleOmega(i, lct(j), true);

    // t1 is the latest time points between est(i) and contact[i]
    auto t1{profile[contact[i]].time < est(i) ? est_shared[i] : contact[i]};

    // available
    auto available{data[lct_shared[j]].available - data[t1].available};
    auto t2{ect_shared[j]};
    auto t3{lct_shared[i]};

    T maxoverflow;
    if (ect(i) < lct(j)) {
      if (available < minenergy(i)) {
        maxoverflow = data[t3].overlap - data[t3].slackUnder -
                      data[t1].overlap + data[t1].slackUnder;
      } else {
        maxoverflow = data[t2].overlap - data[t3].slackUnder -
                      data[t1].overlap + data[t1].slackUnder;
      }
    } else {
      maxoverflow = data[t2].overlap - data[t2].slackUnder - data[t1].overlap +
                    data[t1].slackUnder;
    }

    T adjustment{-Constant::Infinity<T>};
    if (maxoverflow > 0) {
      auto next{profile.at(t1)};
      while (next != profile.end()) {
        auto t{next};
        ++next;
        auto overlap{data[next.index].overlap - data[t.index].overlap};
        if (maxoverflow > overlap) {
          maxoverflow -= overlap;
        } else {
          adjustment = std::min(
              next->time,
              t->time + maxoverflow / (data[t.index].consumption -
                                       capacity.max(solver) + mindemand(i)));
          break;
        }
      }
    }

    computeForwardExplanation(i);
    pruning.push_back(task[i].start.after(adjustment));

    std::cout << " *** adjustment: " << pruning.back() << " because";
    for (auto l : explanation[num_explanations - 1])
      std::cout << " " << l;
    std::cout << std::endl;
  }
}


template <typename T> void CumulativeEdgeFinding<T>::forwardDetection() {

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "\nstart forward detection\n";
  }
#endif

  while (not in_conflict.empty()) {
    prec[in_conflict.back()] = -1;
    in_conflict.pop_back();
  }

  //        auto cap{capacity.max(solver)};
  auto stop{lct_order.rend()};
  --stop;

  // explore the tasks by decreasing lct
  //    for (auto ii{lct_order.rbegin()}; ii != stop;) {
  int k{static_cast<int>(lct_order.size() - 1)};
  while (k >= 0) {
    auto i{lct_order[k]};

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << " - analyse tasks whose lct is " << lct(i) << std::endl;
    }
#endif

    // remove tasks whose lct is larger than or equal to lct(*ii)
    shrinkLeftCutToTime(lct(i));

    // if there are no more tasks, all those in the current level have the same
    // lct and we can stop
    if (leftcut_pointer == 0) {
      break;
    }

    // lct of the task that precedes all the task whose lct is lct(i)
    auto j{lct_order[leftcut_pointer - 1]};

    // tasks in the range [leftcut_pointer, *ii) all have a lct equal to
    // lct(*ii)
    // - add their "prime" versions one by one and run scheduleOmega
    while (k >= leftcut_pointer) {

#ifdef DBG_SEF
      if (DBG_SEF) {
        std::cout << "  * " << prettyTask(i) << ": schedule tasks on [0,"
                  << lct(j) << ")\n";
      }
#endif

      addPrime(i, j);
      scheduleOmega(i, lct(j));
      rmPrime(i, j);

      computeBound(i);

      if (beta != -1) {
#ifdef DBG_SEF
        if (DBG_SEF) {
          std::cout << "  - beta = " << beta << " (task " << task[beta].id()
                    << ")" << std::endl;
        }
#endif

        assert(lct(lct_order[leftcut_pointer]) > lct(beta));
        shrinkLeftCutToTime(lct(beta) + 1);

        addPrime(i, beta);
        auto ect_i_H = scheduleOmega(i, lct(beta));
        rmPrime(i, beta);

        if (ect_i_H > lct(beta)) {

#ifdef DBG_SEF
          if (DBG_SEF) {
            std::cout << " task " << task[i].id()
                      << " is in conflict with task interval ending at time "
                      << lct(beta) << std::endl;
          }
#endif

          prec[i] = beta;
          in_conflict.push_back(i);
        }
      }

      if (prec[i] == -1 and alpha != -1) {

#ifdef DBG_SEF
        if (DBG_SEF) {
          std::cout << "  - alpha = " << alpha << " (task " << task[alpha].id()
                    << ")" << std::endl;
        }
#endif

        assert(lct(lct_order[leftcut_pointer]) > lct(alpha));
        shrinkLeftCutToTime(lct(alpha) + 1);

        //                std::cout << profile << std::endl;

        addPrime(i, alpha);
        auto ect_i_H = scheduleOmega(i, lct(alpha));
        rmPrime(i, alpha);

        if (ect_i_H > lct(alpha)) {

#ifdef DBG_SEF
          if (DBG_SEF) {
            std::cout << " task " << task[i].id()
                      << " is in conflict with task interval "
                      << task[*(lct_order.begin())].id() << ".."
                      << task[alpha].id() << std::endl;
          }
#endif

          prec[i] = alpha;
          in_conflict.push_back(i);
        }
      }

      if (alpha != -1 or beta != -1) {
        growLeftCutToTime(lct(i));
      }

      i = lct_order[--k];
    }
  }
}


template <typename T>
void CumulativeEdgeFinding<T>::computeForwardExplanation(const int i) {
    
    auto h{num_explanations};
    if (explanation.size() <= h) {
      explanation.resize(h + 1);
    } else {
      explanation[h].clear();
    }
    ++num_explanations;
    
    auto k{leftcut_pointer-1};
    auto j{lct_order[k]};

    assert(lct(j) == lct(prec[i]));
    
    auto t{profile[contact[i]].time};
    
    assert(ect(j) > t);
 
    do {
        explanation[h].push_back(task[j].start.after(est(j)));
        explanation[h].push_back(task[j].end.before(lct(j)));
    } while(k > 0 and ect(j = lct_order[--k]) > t);
    explanation[h].push_back(task[i].start.after(est(i)));

}

template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T>, const hint,
                                      std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding"; // "data/sample/j309_5.sm"

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
std::ostream &CumulativeEdgeFinding<T>::print_reason(std::ostream &os, const hint) const {
  os << "cumulative-edge-finding";
  return os;
}

} // namespace tempo

#endif
