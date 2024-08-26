/************************************************
 * Tempo CumulativeEdgeFinding.hpp
 * Implementation of the "strong edge-finding" algorithm as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder rule for the cumulative constraint. In Proceedings of the Twenty-Fifth International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI Press, 3103–3109.
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
    os << "(" << time << "|" << capacity << "|" << increment << "...)";
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

  // helpers
    List<Timepoint<T>> profile;
    std::vector<int> est_;
    std::vector<int> ect_;
    std::vector<int> lct_;
    int sentinel;
    
    std::vector<std::vector<int>> triggers;

    int level;
    SubscriberHandle restartToken;
    SubscriberHandle backtrackToken;
    
    
    void initialiseProfile();

  public:
    static int est_flag;
    static int ect_flag;
    static int lct_flag;
    static int dem_flag;
    
    
  template <typename ItTask, typename ItNVar, typename ItBVar>
  CumulativeEdgeFinding(Solver<T> &solver, const Interval<T> sched, const NumericVar<T> capacity, const ItTask beg_task,
                         const ItTask end_task, const ItNVar beg_dem, const ItBVar beg_disj);
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

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
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
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem,
    const ItBVar beg_disj)
    : m_solver(solver)
 , restartToken(m_solver.SearchRestarted.subscribe_handled(
                            [this](const bool) {
#ifdef DBG_SEF
                                if (DBG_SEF) {
                                    std::cout << "signal restart, re-init profile\n";
                                }
#endif
                                initialiseProfile(); }))
      , backtrackToken(m_solver.BackTrackCompleted.subscribe_handled([this]() {
#ifdef DBG_SEF
          if (DBG_SEF) {
              std::cout << "signal backtrack " << m_solver.level() << "/" << level << "\n";
          }
#endif
          if (level >= m_solver.level()) {
#ifdef DBG_SEF
              if (DBG_SEF) {
                  std::cout << "backtrack over level " << level << ", re-init profile\n";
                  //              level = m_solver.level();
              }
#endif
              initialiseProfile();
          }
      })) 

{

  level = -1;

  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    the_tasks.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  precedence.resize(the_tasks.size());
  est_.resize(the_tasks.size());
  ect_.resize(the_tasks.size());
  lct_.resize(the_tasks.size());

  auto maxcap{capacity.max(m_solver)};
  for (unsigned i = 0; i < the_tasks.size(); ++i) {
    precedence[i].resize(the_tasks.size());
    est_[i] = profile.create_element(est(i), maxcap, mindemand(i), mindemand(i),
                                     0, 0, 0, 0);
    ect_[i] =
        profile.create_element(ect(i), maxcap, -mindemand(i), 0, 0, 0, 0, 0);
    lct_[i] =
        profile.create_element(lct(i), maxcap, 0, -mindemand(i), 0, 0, 0, 0);
  }

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

  std::vector<int> event_ordering(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);

  std::sort(event_ordering.begin(), event_ordering.end(),
            [this](const int i, const int j) {
              return this->profile[i].time < this->profile[j].time;
            });

  auto elt{event_ordering.begin()};
  profile.add_front(*elt);

  //        std::cout << profile << std::endl;

  while (++elt != event_ordering.end()) {
    profile.add_after(*(elt - 1), *elt);

    //            std::cout << "add " <<  *elt << " after " << *(elt-1) <<
    //            std::endl;

    //            std::cout << profile << std::endl;
  }

  sentinel = profile.create_element(profile[*event_ordering.rbegin()].time, 0,
                                    0, 0, 0, 0, 0, 0);

  profile.add_after(*event_ordering.rbegin(), sentinel);

//  std::cout << profile << std::endl;
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}

template <typename T> void CumulativeEdgeFinding<T>::initialiseProfile() {
    
#ifdef DBG_SEF
    if (DBG_SEF) {
        std::cout << "re-init profile @lvl=" << m_solver.level() << std::endl;
    }
#endif
    
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
bool CumulativeEdgeFinding<T>::notify(const Literal<T> 
#ifdef DBG_SEF
                                      l
#endif
                                      , const int r) {
//    int n{static_cast<int>(the_tasks.size())};

#ifdef DBG_SEF
if (DBG_SEF) {
  std::cout << "\nnotify (" << this->id() << ") " << m_solver.pretty(l) << " @lvl="
            << m_solver.level() << "/" << level << std::endl;
}
#endif

//if (level > m_solver.level()) {
//
//#ifdef DBG_SEF
//  if (DBG_SEF) {
//    std::cout << "backtrack, re-initializing\n";
//  }
//#endif
//
//  return true;
//}

#ifdef DBG_SEF

for (unsigned i{0}; i < the_tasks.size(); ++i) {
  if (profile[est_[i]].time > est(i)) {
    std::cout << " (beg not " << this->id() << ") bug profile ahead (est(" << i
              << ")) : " << profile[est_[i]].time << " / " << prettyTask(i)
              << std::endl;
    exit(1);
  }
  if (profile[ect_[i]].time > ect(i)) {
    std::cout << " (beg not " << this->id() << ") bug profile ahead (ect(" << i
              << ")) : " << profile[ect_[i]].time << " / " << prettyTask(i)
              << std::endl;
    exit(1);
  }
  if (profile[lct_[i]].time < lct(i)) {
    std::cout << " (beg not " << this->id() << ") bug profile ahead (lct(" << i
              << ")) : " << profile[lct_[i]].time << " / " << prettyTask(i)
              << std::endl;
    exit(1);
  }
}

#endif
    
    level = m_solver.level();
    
    for(auto t : triggers[r]) {
        auto flag{t%4};
        auto i{t/4};

#ifdef DBG_SEF
        if (DBG_SEF)
          std::cout << "task " << prettyTask(i) << " -- " << est(i) << " "
                    << ect(i) << " " << lct(i) << " -- " << l.value()
                    << std::endl;
#endif
        if(flag == est_flag) {

#ifdef DBG_SEF
          if (DBG_SEF)
            std::cout << "est of task " << the_tasks[i] << std::endl;
#endif

            // est of the_tasks[r] has changed (increased)

#ifdef DBG_SEF
          if (DBG_SEF)
            if (profile[est_[i]].time >= est(i)) {
              std::cout << " bug profile ahead : " << profile[est_[i]].time
                        << std::endl;
            }
#endif

          //          auto t =
          profile[est_[i]].time = est(i);
          //          auto s{profile.at(est_[i])};
          //          auto j{s};
          //          do
          //            ++j;
          //          while (j->time < t);
          //          --j;
          //          if (j != s) {
          //#ifdef DBG_SEF
          //            if (DBG_SEF)
          //              std::cout << " -insert " << profile[est_[i]] << "
          //              after " << *j
          //                        << std::endl;
          ////            std::cout << profile << std::endl;
          //#endif
          //
          //        profile.remove(est_[i]);
          //        profile.add_after(j.index, est_[i]);
          //      }
          //#ifdef DBG_SEF
          //      else if (DBG_SEF)
          //        std::cout << " - " << profile[est_[i]] << "'s rank has not
          //        changed"
          //                  << std::endl;
          //#endif

        } else if (flag == ect_flag) {

#ifdef DBG_SEF
          if (DBG_SEF)
            std::cout << "ect of task " << the_tasks[i] << std::endl;
#endif

#ifdef DBG_SEF
          if (DBG_SEF)
            if (profile[ect_[i]].time >= ect(i)) {
              std::cout << " bug profile ahead : " << profile[ect_[i]].time
                        << std::endl;
            }
#endif

          // ect of the_tasks[r-n] has changed (increased)
          //      auto t =
          profile[ect_[i]].time = ect(i);
          //      auto s{profile.at(ect_[i])};
          //      auto j{s};
          //      do
          //        ++j;
          //      while (j->time < t);
          //      --j;
          //      if (j != s) {
          //#ifdef DBG_SEF
          //        if (DBG_SEF)
          //          std::cout << " -insert " << profile[est_[i]] << " after "
          //          << *j
          //                    << std::endl;
          ////            std::cout << profile << std::endl;
          //#endif
          //
          //        profile.remove(ect_[i]);
          //        profile.add_after(j.index, ect_[i]);
          //      }
          //#ifdef DBG_SEF
          //      else if (DBG_SEF)
          //        std::cout << " - " << profile[ect_[i]] << "'s rank has not
          //        changed"
          //                  << std::endl;
          //#endif

        } else if (flag == lct_flag) {
#ifdef DBG_SEF
          if (DBG_SEF)
            std::cout << "lct of task " << the_tasks[i] << std::endl;
#endif

#ifdef DBG_SEF
          if (DBG_SEF)
            if (profile[lct_[i]].time <= lct(i)) {
              std::cout << " bug profile ahead : " << profile[lct_[i]].time
                        << std::endl;
            }
#endif

          // lct of the_tasks[r-2*n] has changed (decreased)
          //      auto t =
          profile[lct_[i]].time = lct(i);
          //      auto j{profile.at(lct_[i])};
          //
          //      //            std::cout << *j << "\nin profile:\n" << profile
          //      <<
          //      //            std::endl;
          //
          //      do
          //        --j;
          //      while (j->time > t);
          //      if (profile.next(j.index) != lct_[i]) {
          //#ifdef DBG_SEF
          //        if (DBG_SEF)
          //          std::cout << " -insert " << profile[lct_[i]] << " after "
          //          << *j
          //                    << std::endl;
          ////            std::cout << profile << std::endl;
          //#endif
          //
          //        profile.remove(lct_[i]);
          //        profile.add_after(j.index, lct_[i]);
          //      }
          //#ifdef DBG_SEF
          //      else if (DBG_SEF)
          //        std::cout << " - " << profile[lct_[i]] << "'s rank has not
          //        changed"
          //                  << std::endl;
          //#endif
          //
        } else {
          std::cout << "demand of task " << the_tasks[i] << std::endl;
        }
    }

#ifdef DBG_SEF
    if (DBG_SEF) {
      for (auto e{profile.begin()}; e != profile.end(); ++e) {
        std::cout << " " << e.index << ":" << e->time;
      }
      std::cout << std::endl;
    }

    for (unsigned i{0}; i < the_tasks.size(); ++i) {
      if (profile[est_[i]].time > est(i)) {
        std::cout << " (end not) bug profile ahead : " << profile[est_[i]].time
                  << " / " << prettyTask(i) << std::endl;
        exit(1);
      }
      if (profile[ect_[i]].time > ect(i)) {
        std::cout << " (end not) bug profile ahead : " << profile[ect_[i]].time
                  << " / " << prettyTask(i) << std::endl;
        exit(1);
      }
      if (profile[lct_[i]].time < lct(i)) {
        std::cout << " (end not) bug profile ahead : " << profile[lct_[i]].time
                  << " / " << prettyTask(i) << std::endl;
        exit(1);
      }
    }

#endif

    return true;
}


template <typename T> void CumulativeEdgeFinding<T>::propagate() {
    
//    level = m_solver.level();

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "\npropagate " << level << "/" << m_solver.level()
              << std::endl;
    for (auto e{profile.begin()}; e != profile.end(); ++e) {
      std::cout << " " << e.index << ":" << e->time;
    }
    std::cout << std::endl;
    std::cout << "hello " << prettyTask(1) << " " << profile[est_[1]].time
              << ".." << profile[lct_[1]].time << std::endl;
  }
  for (unsigned i{0}; i < the_tasks.size(); ++i) {
    if (profile[est_[i]].time > est(i)) {
      std::cout << " (beg prop) bug profile ahead : " << profile[est_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
    if (profile[ect_[i]].time > ect(i)) {
      std::cout << " (beg prop) bug profile ahead : " << profile[ect_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
    if (profile[lct_[i]].time < lct(i)) {
      std::cout << " (beg prop) bug profile ahead : " << profile[lct_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
  }
#endif

  //    if (level > m_solver.level()) {
  //
  //        for(unsigned i{0}; i<the_tasks.size(); ++i) {
  //            profile[est_[i]].time = est(i);
  //            profile[ect_[i]].time = ect(i);
  //            profile[lct_[i]].time = lct(i);
  //        }
  //
  //#ifdef DBG_SEF
  //        if (DBG_SEF) {
  //            std::cout << "backtrack, re-initializing\n";
  //        }
  //#endif
  //    }

  //    std::cout << profile << std::endl;

  auto e{profile.begin()};
  ++e;
  for (; e != profile.end(); ++e) {

    //        std::cout << "insert " << *e << std::endl;

    auto need_swap{false};
    auto i{e};
    do {
      --i;

      //          std::cout << " - " << *i << "?\n";

      if (i->time <= e->time) {
        break;
      } else {
        need_swap = true;
      }
    } while (i != profile.end());
    if (need_swap) {
      auto x{e};

      //          std::cout << " * swap " << *x << " and " << *i << std::endl;

      profile.remove(x.index);
      profile.add_after(i.index, x.index);
    }
  }

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "reordering\n";
    for (auto e{profile.begin()}; e != profile.end(); ++e) {
      std::cout << " " << e.index << ":" << e->time;
    }
    std::cout << std::endl;
  }
#endif

  //        exit(1);
  //  }
  //  level = m_solver.level();

#ifdef DBG_SEF
  if (DBG_SEF) {

    for (auto e{profile.begin()}; e != profile.end(); ++e) {
      std::cout << " " << e.index << ":" << e->time;
    }
    std::cout << std::endl;
    std::cout << "hello " << prettyTask(1) << " " << profile[est_[1]].time
              << ".." << profile[lct_[1]].time << std::endl;
  }
  T previous{-Constant::Infinity<T>};
  for (auto e{profile.begin()}; e != profile.end(); ++e) {
    if (e->time < previous) {
      std::cout << "not ordered (event " << *e << ") in\n"
                << profile << std::endl;
      exit(1);
    }
    previous = e->time;
  }

  for (unsigned i{0}; i < the_tasks.size(); ++i) {
    if (profile[est_[i]].time > est(i)) {
      std::cout << " (end prop) bug profile ahead : " << profile[est_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
    if (profile[ect_[i]].time > ect(i)) {
      std::cout << " (end prop) bug profile ahead : " << profile[ect_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
    if (profile[lct_[i]].time < lct(i)) {
      std::cout << " (end prop) bug profile ahead : " << profile[lct_[i]].time
                << " / " << prettyTask(i) << std::endl;
      exit(1);
    }
  }
#endif

  //    std::cout << "ok\n";
}



template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T>, const hint,
                                       std::vector<Literal<T>> &) {

 
}

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
