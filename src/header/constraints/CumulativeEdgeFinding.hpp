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
#include <numeric>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Model.hpp"
#include "constraints/Constraint.hpp"
#include "util/List.hpp"

namespace tempo {

template <typename T = int> struct Timepoint {

  Timepoint() {}

  Timepoint(T time, T capacity, T increment, T incrementMax, T overflow,
            T minimumOverflow, T hMaxTotal, T hreal)
      : time(time), capacity(capacity), increment(increment),
        incrementMax(incrementMax), overflow(overflow),
        minimumOverflow(minimumOverflow), hMaxTotal(hMaxTotal), hreal(hreal) {}

  T time{0};
  T capacity{0};
  T increment{0};
  T incrementMax{0};
  T overflow{0};
  T minimumOverflow{0};
  T hMaxTotal{0};
  T hreal{0};

  std::ostream &display(std::ostream &os) const {
    os << "(" << time << "|" << capacity << "|" << increment << "|" << incrementMax << "...)";
    return os;
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Timepoint<T> &x) {
  return x.display(os);
}

template<typename T>
class Solver;

template <typename T> class CumulativeEdgeFinding : public Constraint<T> {
private:
  Solver<T> &m_solver;
  NumericVar<T> capacity;
  Interval<T> schedule;
  std::vector<Interval<T>> the_tasks;
  std::vector<NumericVar<T>> demand;
  std::vector<std::vector<Literal<T>>> precedence;
  std::vector<int> prec;
  std::vector<int> prev;

  // helpers
  List<Timepoint<T>> profile;
  std::vector<int> est_;
  std::vector<int> ect_;
  std::vector<int> lct_;
    
    std::vector<int> event_ordering; //(profile.size());

  std::vector<int> lct_order;
  int sentinel;
    int alpha;
    int beta;

  std::vector<std::vector<int>> triggers;

  //  int level;
  //  SubscriberHandle restartToken;
  //  SubscriberHandle backtrackToken;

  void initialiseProfile();

public:
  static int est_flag;
  static int ect_flag;
  static int lct_flag;
  static int dem_flag;

  template <typename ItTask, typename ItNVar, typename ItBVar>
  CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched,
                        const NumericVar<T> capacity, const ItTask beg_task,
                        const ItTask end_task, const ItNVar beg_dem,
                        const ItBVar beg_disj);
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
  T energy(const unsigned i) const;
    
//    T slack(const unsigned i, const unsigned j, const T e) const;

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  T scheduleOmega();
  T scheduleOmegaMinus(const int b, std::vector<int>::reverse_iterator& j);
  void computeBound(const int i);
  void buildFullProfile();
  void addPrime(const int i);
  void rmPrime();
  void horizontallyElasticEdgeFinder();

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;
};

template <typename T>
int CumulativeEdgeFinding<T>::est_flag = 0;

template <typename T>
int CumulativeEdgeFinding<T>::ect_flag = 1;

template <typename T>
int CumulativeEdgeFinding<T>::lct_flag = 2;

template <typename T>
int CumulativeEdgeFinding<T>::dem_flag = 3;

template <typename T>
std::string CumulativeEdgeFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  ss << the_tasks[i] << ": [" << est(i) << ".." << lct(i) << "] ("
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
  for (auto k{ect(i)}; k < lct(i) - 1; ++k) {
    ss << ".";
  }
  ss << "] " << est(i) << ".." << lct(i);
  return ss.str();
}



template <typename T> T CumulativeEdgeFinding<T>::est(const unsigned i) const {
  return the_tasks[i].getEarliestStart(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lst(const unsigned i) const {
  return the_tasks[i].getLatestStart(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::ect(const unsigned i) const {
  return the_tasks[i].getEarliestEnd(m_solver);
}

template <typename T> T CumulativeEdgeFinding<T>::lct(const unsigned i) const {
  return the_tasks[i].getLatestEnd(m_solver);
}

//template <typename T>
//T CumulativeEdgeFinding<T>::slack(const unsigned i, const unsigned j, const T e) const {
//    // i is the index of a task in 'the_tasks'
//    // j is the index of a task in 'the_tasks' which is the right limit of an LCUT in conflict with i
//    // e is the ect of the profile ending in j
//    if(lct(j) <= ect(j)) {
//        return (capacity.max(m_solver) - demand(i)) * (lct(j) - e);
//    } else {
//        return capacity.max(m_solver) * (lct(j) - e);
//    }
//}

template <typename T>
T CumulativeEdgeFinding<T>::minduration(const unsigned i) const {
  return the_tasks[i].minDuration(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxduration(const unsigned i) const {
  return the_tasks[i].maxDuration(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::mindemand(const unsigned i) const {
  return demand[i].min(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::maxdemand(const unsigned i) const {
  return demand[i].max(m_solver);
}

template <typename T>
T CumulativeEdgeFinding<T>::energy(const unsigned i) const {
  return mindemand(i) * minduration(i);
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem,
    const ItBVar beg_disj)
    : m_solver(solver)
//      ,restartToken(
//          m_solver.SearchRestarted.subscribe_handled([this](const bool) {
//#ifdef DBG_SEF
//            if (DBG_SEF) {
//              std::cout << "signal restart, re-init profile\n";
//            }
//#endif
//            initialiseProfile();
//          }))
//      ,backtrackToken(m_solver.BackTrackCompleted.subscribe_handled([this]() {
//#ifdef DBG_SEF
//        if (DBG_SEF) {
//          std::cout << "signal backtrack " << m_solver.level() << "/" << level
//                    << "\n";
//        }
//#endif
//        if (level >= m_solver.level()) {
//#ifdef DBG_SEF
//          if (DBG_SEF) {
//            std::cout << "backtrack over level " << level
//                      << ", re-init profile\n";
//            //              level = m_solver.level();
//          }
//#endif
//          initialiseProfile();
//        }
//      }))

{

  //  level = -1;

  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    the_tasks.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  auto ip{the_tasks.size()};

  prec.resize(ip);
  precedence.resize(ip);
  est_.resize(ip + 1);
  ect_.resize(ip + 1);
  lct_.resize(ip + 1);

  auto maxcap{capacity.max(m_solver)};

  for (unsigned i = 0; i < ip; ++i) {
    precedence[i].resize(the_tasks.size());
    est_[i] = profile.create_element(est(i), maxcap, mindemand(i), mindemand(i),
                                     0, 0, 0, 0);
    ect_[i] =
        profile.create_element(ect(i), maxcap, -mindemand(i), 0, 0, 0, 0, 0);
    lct_[i] =
        profile.create_element(lct(i), maxcap, 0, -mindemand(i), 0, 0, 0, 0);
  }
  sentinel = profile.create_element(Constant::Infinity<T>, 0,
                                      0, 0, 0, 0, 0, 0);

  auto ep{beg_disj};
  for (auto ip{beg_task}; ip != end_task; ++ip) {
    for (auto jp{ip + 1}; jp != end_task; ++jp) {
      auto x{*ep};

      auto i{std::distance(beg_task, ip)};
      auto j{std::distance(beg_task, jp)};
      precedence[i][j] = m_solver.boolean.getLiteral(false, x);

      x = *(++ep);

      precedence[j][i] = m_solver.boolean.getLiteral(false, x);
      ++ep;
    }
  }

  lct_order.resize(the_tasks.size());
  std::iota(lct_order.begin(), lct_order.end(), 0);
  //          std::sort(lct_order.begin(), lct_order.end(),
  //                    [this](const int i, const int j) {
  //                      return lct(i) < lct(j);
  //                    });

//  std::vector<int> event_ordering(profile.size());

  prev.resize(profile.size());
  event_ordering.resize(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);
    est_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0);
    ect_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0);
    lct_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0);
//  std::sort(event_ordering.begin(), event_ordering.end(),
//            [this](const int i, const int j) {
//              return this->profile[i].time < this->profile[j].time;
//            });

  //  auto elt{event_ordering.begin()};
  //  profile.add_front(*elt);

  //  //        std::cout << profile << std::endl;
  //
  //  while (++elt != event_ordering.end()) {
  //    profile.add_after(*(elt - 1), *elt);
  //
  //    //            std::cout << "add " <<  *elt << " after " << *(elt-1) <<
  //    //            std::endl;
  //
  //    //            std::cout << profile << std::endl;
  //  }

  

  //  profile.add_after(*event_ordering.rbegin(), sentinel);

  //  std::cout << profile << std::endl;
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}

template <typename T> void CumulativeEdgeFinding<T>::initialiseProfile() {

  //#ifdef DBG_SEF
  //  if (DBG_SEF) {
  //    std::cout << "re-init profile @lvl=" << m_solver.level() << std::endl;
  //  }
  //#endif

  for (unsigned i{0}; i < the_tasks.size(); ++i) {
    profile[est_[i]].time = est(i);
    profile[ect_[i]].time = ect(i);
    profile[lct_[i]].time = lct(i);
  }
}

template <typename T> void CumulativeEdgeFinding<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < the_tasks.size(); ++i) {
    auto k{m_solver.wake_me_on(lb<T>(the_tasks[i].start.id()), this->id())};
    if (triggers.size() <= k)
      triggers.resize(k + 1);
    triggers[k].push_back(est_flag + 4 * i);
  }
  for (size_t i{0}; i < the_tasks.size(); ++i) {
    auto k{m_solver.wake_me_on(lb<T>(the_tasks[i].end.id()), this->id())};
    if (triggers.size() <= k)
      triggers.resize(k + 1);
    triggers[k].push_back(ect_flag + 4 * i);
  }
  for (size_t i{0}; i < the_tasks.size(); ++i) {
    auto k{m_solver.wake_me_on(ub<T>(the_tasks[i].end.id()), this->id())};
    if (triggers.size() <= k)
      triggers.resize(k + 1);
    triggers[k].push_back(lct_flag + 4 * i);
  }
  for (size_t i{0}; i < the_tasks.size(); ++i) {
    auto k{m_solver.wake_me_on(lb<T>(demand[i].id()), this->id())};
    if (triggers.size() <= k)
      triggers.resize(k + 1);
    triggers[k].push_back(dem_flag + 4 * i);
  }
}

template <typename T>
bool CumulativeEdgeFinding<T>::notify(const Literal<T>, const int) {
  return true;
}

//
//#ifdef DBG_SEF
// if (DBG_SEF) {
//  std::cout << "\nnotify (" << this->id() << ") " << m_solver.pretty(l)
//            << " @lvl=" << m_solver.level() << "/" << level << std::endl;
//}
//#endif
//
//
//#ifdef DBG_SEF
//
// for (unsigned i{0}; i < the_tasks.size(); ++i) {
//  if (profile[est_[i]].time > est(i)) {
//    std::cout << " (beg not " << this->id() << ") bug profile ahead (est(" <<
//    i
//              << ")) : " << profile[est_[i]].time << " / " << prettyTask(i)
//              << std::endl;
//    exit(1);
//  }
//  if (profile[ect_[i]].time > ect(i)) {
//    std::cout << " (beg not " << this->id() << ") bug profile ahead (ect(" <<
//    i
//              << ")) : " << profile[ect_[i]].time << " / " << prettyTask(i)
//              << std::endl;
//    exit(1);
//  }
//  if (profile[lct_[i]].time < lct(i)) {
//    std::cout << " (beg not " << this->id() << ") bug profile ahead (lct(" <<
//    i
//              << ")) : " << profile[lct_[i]].time << " / " << prettyTask(i)
//              << std::endl;
//    exit(1);
//  }
//}
//
//#endif
//
// level = m_solver.level();
//
// for (auto t : triggers[r]) {
//  auto flag{t % 4};
//  auto i{t / 4};
//
//#ifdef DBG_SEF
//  if (DBG_SEF)
//    std::cout << "task " << prettyTask(i) << " -- " << est(i) << " " << ect(i)
//              << " " << lct(i) << " -- " << l.value() << std::endl;
//#endif
//  if (flag == est_flag) {
//
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      std::cout << "est of task " << the_tasks[i] << std::endl;
//#endif
//
//      // est of the_tasks[r] has changed (increased)
//
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      if (profile[est_[i]].time >= est(i)) {
//        std::cout << " bug profile ahead : " << profile[est_[i]].time
//                  << std::endl;
//      }
//#endif
//
//    //          auto t =
//    profile[est_[i]].time = est(i);
//    //          auto s{profile.at(est_[i])};
//    //          auto j{s};
//    //          do
//    //            ++j;
//    //          while (j->time < t);
//    //          --j;
//    //          if (j != s) {
//    //#ifdef DBG_SEF
//    //            if (DBG_SEF)
//    //              std::cout << " -insert " << profile[est_[i]] << "
//    //              after " << *j
//    //                        << std::endl;
//    ////            std::cout << profile << std::endl;
//    //#endif
//    //
//    //        profile.remove(est_[i]);
//    //        profile.add_after(j.index, est_[i]);
//    //      }
//    //#ifdef DBG_SEF
//    //      else if (DBG_SEF)
//    //        std::cout << " - " << profile[est_[i]] << "'s rank has not
//    //        changed"
//    //                  << std::endl;
//    //#endif
//
//  } else if (flag == ect_flag) {
//
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      std::cout << "ect of task " << the_tasks[i] << std::endl;
//#endif
//
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      if (profile[ect_[i]].time >= ect(i)) {
//        std::cout << " bug profile ahead : " << profile[ect_[i]].time
//                  << std::endl;
//      }
//#endif
//
//    // ect of the_tasks[r-n] has changed (increased)
//    //      auto t =
//    profile[ect_[i]].time = ect(i);
//    //      auto s{profile.at(ect_[i])};
//    //      auto j{s};
//    //      do
//    //        ++j;
//    //      while (j->time < t);
//    //      --j;
//    //      if (j != s) {
//    //#ifdef DBG_SEF
//    //        if (DBG_SEF)
//    //          std::cout << " -insert " << profile[est_[i]] << " after "
//    //          << *j
//    //                    << std::endl;
//    ////            std::cout << profile << std::endl;
//    //#endif
//    //
//    //        profile.remove(ect_[i]);
//    //        profile.add_after(j.index, ect_[i]);
//    //      }
//    //#ifdef DBG_SEF
//    //      else if (DBG_SEF)
//    //        std::cout << " - " << profile[ect_[i]] << "'s rank has not
//    //        changed"
//    //                  << std::endl;
//    //#endif
//
//  } else if (flag == lct_flag) {
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      std::cout << "lct of task " << the_tasks[i] << std::endl;
//#endif
//
//#ifdef DBG_SEF
//    if (DBG_SEF)
//      if (profile[lct_[i]].time <= lct(i)) {
//        std::cout << " bug profile ahead : " << profile[lct_[i]].time
//                  << std::endl;
//      }
//#endif
//
//    // lct of the_tasks[r-2*n] has changed (decreased)
//    //      auto t =
//    profile[lct_[i]].time = lct(i);
//    //      auto j{profile.at(lct_[i])};
//    //
//    //      //            std::cout << *j << "\nin profile:\n" << profile
//    //      <<
//    //      //            std::endl;
//    //
//    //      do
//    //        --j;
//    //      while (j->time > t);
//    //      if (profile.next(j.index) != lct_[i]) {
//    //#ifdef DBG_SEF
//    //        if (DBG_SEF)
//    //          std::cout << " -insert " << profile[lct_[i]] << " after "
//    //          << *j
//    //                    << std::endl;
//    ////            std::cout << profile << std::endl;
//    //#endif
//    //
//    //        profile.remove(lct_[i]);
//    //        profile.add_after(j.index, lct_[i]);
//    //      }
//    //#ifdef DBG_SEF
//    //      else if (DBG_SEF)
//    //        std::cout << " - " << profile[lct_[i]] << "'s rank has not
//    //        changed"
//    //                  << std::endl;
//    //#endif
//    //
//  } else {
//    std::cout << "demand of task " << the_tasks[i] << std::endl;
//  }
//}
//
//#ifdef DBG_SEF
// if (DBG_SEF) {
//  for (auto e{profile.begin()}; e != profile.end(); ++e) {
//    std::cout << " " << e.index << ":" << e->time;
//  }
//  std::cout << std::endl;
//}
//
// for (unsigned i{0}; i < the_tasks.size(); ++i) {
//  if (profile[est_[i]].time > est(i)) {
//    std::cout << " (end not) bug profile ahead : " << profile[est_[i]].time
//              << " / " << prettyTask(i) << std::endl;
//    exit(1);
//  }
//  if (profile[ect_[i]].time > ect(i)) {
//    std::cout << " (end not) bug profile ahead : " << profile[ect_[i]].time
//              << " / " << prettyTask(i) << std::endl;
//    exit(1);
//  }
//  if (profile[lct_[i]].time < lct(i)) {
//    std::cout << " (end not) bug profile ahead : " << profile[lct_[i]].time
//              << " / " << prettyTask(i) << std::endl;
//    exit(1);
//  }
//}
//
//#endif
//
// return true;
//}

template <typename T> void CumulativeEdgeFinding<T>::buildFullProfile() {
    
//    std::sort(lct_order.begin(), lct_order.end(),
//              [this](const int i, const int j) { return lct(i) < lct(j); });
    
      std::sort(event_ordering.begin(), event_ordering.end(),
                [this](const int i, const int j) {
                  return this->profile[i].time < this->profile[j].time;
                });
    
//    std::cout << profile << std::endl;
    
    auto prev{List<Timepoint<T>>::tail};
    for(auto e : event_ordering) {
        profile.add_after(prev, e);
        prev = e;
        
//        std::cout << " add (" << e << ") " << profile[e] << std::endl;
    }
    profile.add_after(prev, sentinel);
//    
//    std::cout << profile << std::endl;
//    exit(1);
    
}

template <typename T> void CumulativeEdgeFinding<T>::propagate() {

  std::sort(lct_order.begin(), lct_order.end(),
            [this](const int i, const int j) { return lct(i) < lct(j); });

//  profile.add_front(sentinel);
    
    
#ifdef DBG_SEF
        if (DBG_SEF) {
            std::cout << "\n\nstart propagation\n";
                for (auto j : lct_order) {
                    std::cout << "task " << j << ": " << asciiArt(j) << std::endl;
                }
        }
#endif
    

//  int k{0};
    for (auto i : lct_order) {
        profile[est_[i]].time = est(i);
        profile[ect_[i]].time = ect(i);
        profile[lct_[i]].time = lct(i);
//        ++k;
    }
    
    horizontallyElasticEdgeFinder();
//    buildFullProfile();

//    k = 0;
//  for (auto i : lct_order) {
//    auto e{est_[i]};
//    auto j{profile.begin()};
//    auto prev{j};
//    --prev;
//    do {
//      if (j == profile.end() or j->time > profile[e].time) {
//        profile.add_after(prev.index, e);
//        j = profile.at(e);
//        if (e == ect_[i])
//          break;
//        e = ect_[i];
//      } else {
//        prev = j;
//        ++j;
//      }
//    } while (true);
//    profile.add_before(sentinel, lct_[i]);
//    if (k > 0) {
//      auto ect_omega{scheduleOmega()};
//      auto j{lct_order[k - 1]};
//      if (ect_omega > lct(j)) {
//        std::cout << "pruning (" << ect_omega << " > " << lct(j) << ")\n";
//      }
//      for (auto p : profile) {
//        p.capacity = capacity.max(m_solver);
//      }
//    }
//    ++k;
//  }


  profile.clear();
}


template <typename T>
void CumulativeEdgeFinding<T>::addPrime(const int i) {
    auto _est{std::max(est(i), profile.begin()->time)};
    auto _lct{std::min(profile.rbegin()->time, ect(i))};
//    auto _dur{lct-est};
    
    auto maxcap{capacity.max(m_solver)};
    auto ip{the_tasks.size()};
    
//      est_[ip] = profile.create_element(_est, maxcap, mindemand(i), mindemand(i), 0, 0, 0, 0);
//      ect_[ip] = profile.create_element(_lct, maxcap, -mindemand(i), 0, 0, 0, 0, 0);
//      lct_[ip] = profile.create_element(_lct, maxcap, 0, -mindemand(i), 0, 0, 0, 0);
    
    profile[est_[ip]].capacity = profile[ect_[ip]].capacity = profile[lct_[ip]].capacity = maxcap;
    profile[est_[ip]].time = _est;
    profile[ect_[ip]].time = profile[lct_[ip]].time = _lct;
    profile[est_[ip]].increment = profile[est_[ip]].incrementMax = mindemand(i);
    profile[ect_[ip]].increment = profile[lct_[ip]].incrementMax = -mindemand(i);
    
//    std::cout << "create " << est_.back() << "|" << ect_.back() << "|" << lct_.back() << std::endl;
    
    auto p{profile.rbegin()};
    while(p!=profile.rend()) {
        
//        std::cout << " -- " << p->time << "/" << _lct << std::endl;
        
        if(p->time <= _lct) {
            profile.add_after(p.index, lct_[ip]);
            profile.add_after(p.index, ect_[ip]);
            break;
        }
        ++p;
        
//        std::cout << " next=" << p.index << ":" << p->time << std::endl;
    }
    while(p!=profile.rend()) {
        if(p->time <= _est) {
            profile.add_after(p.index, est_[ip]);
            break;
        }
        ++p;
    }
    if(p==profile.rend()) {
        profile.add_front(est_[ip]);
    }
    
//    std::cout << profile << std::endl;
//    exit(1);
}

template <typename T>
void CumulativeEdgeFinding<T>::rmPrime() {
//    auto ip{the_tasks.size()};
    
//    std::cout << "rm " << est_.back() << "|" << ect_.back() << "|" << lct_.back() << std::endl;
//    std::cout << profile << std::endl;
    
    profile.remove(est_.back());
    profile.remove(ect_.back());
    profile.remove(lct_.back());
    
    
    
}


template <typename T>
void CumulativeEdgeFinding<T>::horizontallyElasticEdgeFinder() {
    
    buildFullProfile();
    
    auto cap{capacity.max(m_solver)};
    
    auto stop{lct_order.rend()};
    --stop;
    
    // explore the tasks by decreasing lct
    for(auto ii{lct_order.rbegin()}; ii!=stop;) {
        auto is{ii};
        
        
#ifdef DBG_SEF
        if (DBG_SEF) {
            std::cout << std::endl << profile << " - analyse tasks whose lct is " << lct(*ii) << std::endl; //<< " (remove and add the 'prime'):\n";
        }
#endif
        
        // remove all tasks whose lct is equal the current max
        do {
            
#ifdef DBG_SEF
          if (DBG_SEF) {
            std::cout << "  * rm " << (*ii) << std::endl;
          }
#endif
            
            profile.remove(lct_[*ii]);
            profile.remove(ect_[*ii]);
            profile.remove(est_[*ii]);
            ++ii;

        } while(ii != lct_order.rend() and lct(*is) == lct(*ii));

        // if there are no more tasks, all those in the current level have the same lct and we can stop
        if(ii == lct_order.rend())
            break;
        
        // lct of the task that precedes all the task whose lct is lct(i)
        auto j{*ii};
        auto lct_j{lct(j)};

        // otherwise, add their "prime" versions one by one and run scheduleOmega
        while(is != ii) {
            
            auto i{*is};

            if (ect(i) < lct(i)) {

#ifdef DBG_SEF
              if (DBG_SEF) {
                std::cout << "  * add " << i << "'\n";
              }
#endif
              addPrime(i);
#ifdef DBG_SEF
              if (DBG_SEF) {
                std::cout << profile; //<< std::endl;
              }
#endif

              auto omega_ect{scheduleOmega()};
#ifdef DBG_SEF
              if (DBG_SEF) {
                std::cout << " ect^H = " << omega_ect << " / lct(S) = " << lct_j
                          << std::endl;
                //                    std::cout << "  * rm " << i << "'\n";
              }
#endif
              //                rmPrime();
              if (omega_ect > lct_j) {

#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout << " task " << i
                            << " is in conflict with task interval 0.." << j
                            << std::endl;
                }
#endif
                prec[i] = j;
              } else {

#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout << " compute bound\n";
                }
#endif

                  auto ti{ii};
//                int b = -1;
                computeBound(i);
                if (beta != -1) {
#ifdef DBG_SEF
                  if (DBG_SEF) {
                    std::cout << "  - beta = " << beta << std::endl;
                  }
#endif

//                  b = beta;
                  auto ect_i_H{scheduleOmegaMinus(beta, ti)};

                  if (ect_i_H > lct(beta)) {

#ifdef DBG_SEF
                    if (DBG_SEF) {
                      std::cout << " task " << i
                                << " is in conflict with task interval 0.."
                                << beta << std::endl;
                    }
#endif
                    prec[i] = beta;
                  }
                }
                if (prec[i] == -1 and alpha != -1) {

#ifdef DBG_SEF
                  if (DBG_SEF) {
                    std::cout << "  - alpha = " << alpha << std::endl;
                  }
#endif

//                  b = alpha;
                  auto ect_i_H{
                      scheduleOmegaMinus(alpha, ti)};

                  if (ect_i_H > lct(alpha)) {

#ifdef DBG_SEF
                    if (DBG_SEF) {
                      std::cout << " task " << i
                                << " is in conflict with task interval 0.."
                                << alpha << std::endl;
                    }
#endif
                  }
                }

                  for (; ti != ii;) {

                    --ti;
                    std::cout << "readd " << *ti << std::endl;

                    profile.add_after(prev[est_[*ti]], est_[*ti]);
                    profile.add_after(prev[ect_[*ti]], ect_[*ti]);
                    profile.add_after(prev[lct_[*ti]], lct_[*ti]);
                  }
                }

#ifdef DBG_SEF
              if (DBG_SEF) {
                std::cout << "  * rm " << i << "'\n";
              }
#endif
              rmPrime();

              for (auto p{profile.begin()}; p != profile.end(); ++p) {

                std::cout << *p << std::endl;

                p->capacity = cap;
              }
            }

            ++is;
        }
    }
}


template <typename T>
void CumulativeEdgeFinding<T>::computeBound(const int i) {
    T E{0};
    alpha = -1;
    beta = -1;
    T minSlack[2] = {Constant::Infinity<T>, Constant::Infinity<T>};
    for(auto j : lct_order) {
        if(lct(j) == lct(i))
            break;
        E += energy(j);
        if(lct(j) <= ect(i) and est(i) < lct(j)) {
          auto slack{(capacity.max(m_solver) - mindemand(i)) * lct(j) - E};
          if (slack < minSlack[0] and profile[lct_[j]].overflow > 0) {
            minSlack[0] = slack;
            alpha = j;
          }
        } else if(lct(j) > ect(i)) {
            auto slack{capacity.max(m_solver) * lct(j) - E};
            if(slack < minSlack[1] and profile[lct_[j]].overflow > 0) {
                minSlack[1] = slack;
                beta = j;
            }
        }
    }
}

// template <typename Iter>

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmegaMinus(
    const int b, std::vector<int>::reverse_iterator& jj) {

  //    std::cout << "\nhere:\n" << profile << std::endl;

//  auto jj{j};
  for (; *jj != b; ++jj) {

    //      std::cout << "remove " << *jj << std::endl;

    prev[lct_[*jj]] = profile.prev(lct_[*jj]);
    profile.remove(lct_[*jj]);

    prev[ect_[*jj]] = profile.prev(ect_[*jj]);
    profile.remove(ect_[*jj]);

    prev[est_[*jj]] = profile.prev(est_[*jj]);
    profile.remove(est_[*jj]);
  }

  //    addPrime(i);

  //    std::cout << profile << std::endl;

  auto ect_i_H{scheduleOmega()};

  //    rmPrime();
  //
  //  for (; jj != j;) {
  //
  //      --jj;
  ////      std::cout << "readd " << *jj << std::endl;
  //
  //    profile.add_after(prev[est_[*jj]], est_[*jj]);
  //    profile.add_after(prev[ect_[*jj]], ect_[*jj]);
  //    profile.add_after(prev[lct_[*jj]], lct_[*jj]);
  //  }

  return ect_i_H;
}

template <typename T> T CumulativeEdgeFinding<T>::scheduleOmega() {
    
    
#ifdef DBG_SEF
        if (DBG_SEF) {
            std::cout << "[schedule tasks]" << std::endl;
        }
#endif
    

  auto saved_size{profile.size()};
  auto sentinel{profile.end()};
  --sentinel;

  auto C{capacity.max(m_solver)};
  auto next{profile.begin()};
  T overflow{0};
  T omega_ect{-Constant::Infinity<T>};
  T S{0};
  T h_req{0};

  while (next != sentinel) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
        std::cout << "jump to t_" << t->time << ":";
    }
#endif

    t->overflow = overflow;
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

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << " h_max=" << h_max << ", h_req=" << h_req
                << ", h_cons=" << h_cons
         << ", ov=" << overflow ;
        if(overflow > 0)
        {
            std::cout << "->" << overflow - ((h_cons - h_req) * l) << " @t=" << next->time;
        }
    }
#endif

    // there is some overflow, and it will be resorbed by the next time point
    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
      // then we create a new time point for the moment it will be resorbed
      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
      // h_req)));
      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
      auto new_event{
          profile.create_element(t->time + l, t->capacity, 0, 0, 0, 0, 0, 0)};
      profile.add_after(t.index, new_event);
      next = profile.at(new_event);
    }
    // overflow is the deficit on resource for that period (because tasks are
    // set to their earliest)profile[est_[ip]].time = _est;
    overflow += (h_req - h_cons) * l;
    // once there
    t->capacity = C - h_cons;
    if (overflow > 0)
      omega_ect = Constant::Infinity<T>;
    else if (t->capacity < C)
      omega_ect = profile[profile.next(t.index)].time;

#ifdef DBG_SEF
    if (DBG_SEF) {
      if (omega_ect != -Constant::Infinity<T>)
        std::cout << ", ect=" << omega_ect;
      std::cout << std::endl;
    }
#endif
  }

  while (profile.size() > saved_size) {
    profile.pop_back();
    //        profile.remove(profile.size());
  }
    
#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "return ect^H = " << omega_ect << std::endl;
    }
#endif

  return omega_ect;
}

//
//  //    level = m_solver.level();
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "\npropagate " << level << "/" << m_solver.level()
//              << std::endl;
//    for (auto e{profile.begin()}; e != profile.end(); ++e) {
//      std::cout << " " << e.index << ":" << e->time;
//    }
//    std::cout << std::endl;
//  }
//  for (unsigned i{0}; i < the_tasks.size(); ++i) {
//    if (profile[est_[i]].time > est(i)) {
//      std::cout << " (beg prop) bug profile ahead : " << profile[est_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//    if (profile[ect_[i]].time > ect(i)) {
//      std::cout << " (beg prop) bug profile ahead : " << profile[ect_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//    if (profile[lct_[i]].time < lct(i)) {
//      std::cout << " (beg prop) bug profile ahead : " << profile[lct_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//  }
//#endif
//
//  //    if (level > m_solver.level()) {
//  //
//  //        for(unsigned i{0}; i<the_tasks.size(); ++i) {
//  //            profile[est_[i]].time = est(i);
//  //            profile[ect_[i]].time = ect(i);
//  //            profile[lct_[i]].time = lct(i);
//  //        }
//  //
//  //#ifdef DBG_SEF
//  //        if (DBG_SEF) {
//  //            std::cout << "backtrack, re-initializing\n";
//  //        }
//  //#endif
//  //    }
//
//  //    std::cout << profile << std::endl;
//
//  auto e{profile.begin()};
//  ++e;
//  for (; e != profile.end(); ++e) {
//
//    //        std::cout << "insert " << *e << std::endl;
//
//    auto need_swap{false};
//    auto i{e};
//    do {
//      --i;
//
//      //          std::cout << " - " << *i << "?\n";
//
//      if (i->time <= e->time) {
//        break;
//      } else {
//        need_swap = true;
//      }
//    } while (i != profile.end());
//    if (need_swap) {
//      auto x{e};
//
//      //          std::cout << " * swap " << *x << " and " << *i << std::endl;
//
//      profile.remove(x.index);
//      profile.add_after(i.index, x.index);
//    }
//  }
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "reordering\n";
//    for (auto e{profile.begin()}; e != profile.end(); ++e) {
//      std::cout << " " << e.index << ":" << e->time;
//    }
//    std::cout << std::endl;
//  }
//#endif
//
//  //        exit(1);
//  //  }
//  //  level = m_solver.level();
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//
//    for (auto e{profile.begin()}; e != profile.end(); ++e) {
//      std::cout << " " << e.index << ":" << e->time;
//    }
//    std::cout << std::endl;
//  }
//  T previous{-Constant::Infinity<T>};
//  for (auto e{profile.begin()}; e != profile.end(); ++e) {
//    if (e->time < previous) {
//      std::cout << "not ordered (event " << *e << ") in\n"
//                << profile << std::endl;
//      exit(1);
//    }
//    previous = e->time;
//  }
//
//  for (unsigned i{0}; i < the_tasks.size(); ++i) {
//    if (profile[est_[i]].time > est(i)) {
//      std::cout << " (end prop) bug profile ahead : " << profile[est_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//    if (profile[ect_[i]].time > ect(i)) {
//      std::cout << " (end prop) bug profile ahead : " << profile[ect_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//    if (profile[lct_[i]].time < lct(i)) {
//      std::cout << " (end prop) bug profile ahead : " << profile[lct_[i]].time
//                << " / " << prettyTask(i) << std::endl;
//      exit(1);
//    }
//  }
//
//  std::vector<int> lct_order(the_tasks.size());
//  std::iota(lct_order.begin(), lct_order.end(), 0);
//  std::sort(lct_order.begin(), lct_order.end(),
//            [&](const int x, const int y) { return lct(x) < lct(y); });
//  for (auto i : lct_order) {
//    std::cout << asciiArt(i) << std::endl;
//  }
//
//#endif
//
//  auto saved_size{profile.size()};
//  auto sentinel{profile.end()};
//  --sentinel;
//
//  auto C{capacity.max(m_solver)};
//  auto next{profile.begin()};
//  T overflow{0};
//  T omega_ect{-Constant::Infinity<T>};
//  T S{0};
//  T h_req{0};
//
//  while (next != sentinel) {
//    auto t{next};
//    ++next;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << "jump to t_" << t->time << ": ov=" << overflow;
//    }
//#endif
//
//    t->overflow = overflow;
//    auto l = next->time - t->time;
//
//    // S is the sum of demands of the tasks that could be processed at time
//    S += t->incrementMax;
//    // h_max is the min between the resource's capacity and the total of the
//    // demands
//    auto h_max{std::min(S, C)};
//    // h_req is the total demand counting tasks processed at their earliest
//    h_req += t->increment;
//    // h_cons is the amount of resource actually used in the optimistic
//    scenario
//    // (min between what is required + due from earlier, and what is
//    available) auto h_cons{std::min(h_req + overflow, h_max)};
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << ", h_max=" << h_max << ", h_req=" << h_req
//                << ", h_cons=" << h_cons;
//    }
//#endif
//
//    // there is some overflow, and it will be resorbed by the next time point
//    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
//      // then we create a new time point for the moment it will be resorbed
//      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
//      auto new_event{
//          profile.create_element(t->time + l, t->capacity, 0, 0, 0, 0, 0, 0)};
//      profile.add_after(t.index, new_event);
//    }
//    // overflow is the deficit on resource for that period (because tasks are
//    // set to their earliest)
//    overflow += (h_req - h_cons) * l;
//    // once there
//    t->capacity = C - h_cons;
//    if (t->capacity < C)
//      omega_ect = profile[profile.next(t.index)].time;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      if (omega_ect != -Constant::Infinity<T>)
//        std::cout << ", ect=" << omega_ect;
//      std::cout << std::endl;
//    }
//#endif
//  }
//
//  while (profile.size() > saved_size) {
//    profile.pop_back();
//    //        profile.remove(profile.size());
//  }
//}

template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T>, const hint,
                                      std::vector<Literal<T>> &) {}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding";

#ifdef DBG_EDGEFINDING
  os << "[" << this->id() << "]";
#endif

  os << "(";
  for (auto &t : the_tasks) {
    std::cout << " t" << t.id();
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::print_reason(std::ostream &os,
                                                     const hint) const {
  os << "cumulative-edge-finding";
  return os;
}

} // namespace tempo

#endif
