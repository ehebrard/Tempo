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

  //  Timepoint(T time, T capacity, T increment, T incrementMax, T overflow,
  //            T slackUnder, T available, T overlap)
  //      : time(time), capacity(capacity), increment(increment),
  //        incrementMax(incrementMax), overflow(overflow),
  //    slackUnder(slackUnder), available(available), overlap(overlap) {}

  Timepoint(T time, T capacity, T increment, T incrementMax)
      : time(time), capacity(capacity), increment(increment),
        incrementMax(incrementMax) {}

  T time{0};
  T capacity{0};
  T increment{0};
  T incrementMax{0};
  T overflow{0};
  T slackUnder{0};
  T available{0};
  T overlap{0};
  bool contact{false};

  std::ostream &display(std::ostream &os) const {
    os << "(" << time << "|" << capacity << "|" << increment << "|"
       << incrementMax << "...)";
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
  std::vector<int> prec;   // mapping from the tasks that need to be adjusted to
                           // the corresponding left-cut interval
  SparseSet<> in_conflict; // tasks that have been found to be in conflict with
                           // an left-cut interval
  std::vector<bool> inLeftCut;
  std::vector<int>::reverse_iterator
      lc_ptr; // right side of the left-cut interval

  // helpers
  List<Timepoint<T>> profile;
  std::vector<int> est_;
  std::vector<int> ect_;
  std::vector<int> lct_;

  std::vector<int> event_ordering;

  std::vector<int> lct_order;
  //  int sentinel;
  int alpha;
  int beta;

  T overflow;

  std::vector<std::vector<int>> triggers;

  //  int level;
  //  SubscriberHandle restartToken;
  //  SubscriberHandle backtrackToken;

  //  void initialiseProfile();

public:
  static int est_flag;
  static int ect_flag;
  static int lct_flag;
  static int dem_flag;

  // for (auto t : triggers[r]) {
  //  auto flag{t % 4};
  //  auto i{t / 4};
  // }

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

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  T scheduleOmega(const T C, const T max_lct);
  T scheduleOmegaAdjustment(int u, const T max_lct);
  void computeBound(const int i);
  void buildFullProfile();

  void addTask(const int i);
  void rmTask(const int i);

  void leftCut(const int i);

  void addPrime(const int i);
  void rmPrime();

  void horizontallyElasticEdgeFinderForward();
  void forwardDetection();
  void forwardAdjustment();

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;

    std::vector<bool> inprof; //(profile.size(), false);
    bool verify() {
        
//        std::cout << "hello\n";
        
//        std::vector<bool> inprof(profile.size(), false);
        inprof.clear();
        inprof.resize(profile.size(), false);
        auto prev_t{-Constant::Infinity<T>};
        for(auto p{profile.begin()}; p!=profile.end(); ++p) {
            if(inprof[p.index])
                return false;
            
            inprof[p.index] = true;
            
            if(prev_t > p->time)
                return false;
        }
        
        return checkLeftCut();
    }
    
  bool checkLeftCut() {
//    std::vector<bool> inprof(profile.size(), false);
      int count{0};
    for (auto p{profile.begin()}; p != profile.end(); ++p) {
      if (/*p.index != sentinel or*/ p.index != est_.back() or
          p.index != ect_.back() or p.index != lct_.back())
        ++count;
//      inprof[p.index] = true;
    }
      bool ok{true};
    for (unsigned i{0}; ok and i < the_tasks.size(); ++i) {
      if (inLeftCut[i]) {
        count -= 3;
        if (not(inprof[est_[i]] and inprof[ect_[i]] and inprof[lct_[i]]))
          ok = false;
      }
    }
    ok ^= (count == 0);
    return ok;
  }
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
    if(lct(i) == Constant::Infinity<T>) {
        ss << "=... " << est(i) << "...";
    } else {
        for (auto k{ect(i)}; k < lct(i) - 1; ++k) {
            ss << ".";
        }
        ss << "] " << est(i) << ".." << lct(i);
    }
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
T CumulativeEdgeFinding<T>::energy(const unsigned i) const {
  return mindemand(i) * minduration(i);
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem,
    const ItBVar beg_disj)
    : m_solver(solver) {
  schedule = sched, capacity = cap;

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    the_tasks.push_back(*jp);
    demand.push_back(*dp);
    ++dp;
  }

  auto ip{the_tasks.size()};
  inLeftCut.resize(ip, false);

  in_conflict.reserve(ip);
  prec.resize(ip);
  precedence.resize(ip);
  est_.resize(ip + 1);
  ect_.resize(ip + 1);
  lct_.resize(ip + 1);

  auto maxcap{capacity.max(m_solver)};

  for (unsigned i = 0; i < ip; ++i) {
    precedence[i].resize(the_tasks.size());
    est_[i] =
        profile.create_element(est(i), maxcap, mindemand(i), mindemand(i));
    ect_[i] = profile.create_element(ect(i), maxcap, -mindemand(i), 0);
    lct_[i] = profile.create_element(lct(i), maxcap, 0, -mindemand(i));
  }
  //  sentinel = profile.create_element(Constant::Infinity<T>, 0, 0, 0, 0, 0, 0,
  //  0);

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

  event_ordering.resize(profile.size());
  std::iota(event_ordering.begin(), event_ordering.end(), 1);
  est_[ip] = profile.create_element(0, 0, 0, 0);
  ect_[ip] = profile.create_element(0, 0, 0, 0);
  lct_[ip] = profile.create_element(0, 0, 0, 0);
}

template <typename T> CumulativeEdgeFinding<T>::~CumulativeEdgeFinding() {}

// template <typename T> void CumulativeEdgeFinding<T>::initialiseProfile() {
//   for (unsigned i{0}; i < the_tasks.size(); ++i) {
//     profile[est_[i]].time = est(i);
//     profile[ect_[i]].time = ect(i);
//     profile[lct_[i]].time = lct(i);
//   }
// }

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

template <typename T> void CumulativeEdgeFinding<T>::buildFullProfile() {


    profile.clear();
    
  for (auto i : lct_order) {
    profile[est_[i]].time = est(i);
    profile[ect_[i]].time = ect(i);
    profile[lct_[i]].time = lct(i);
  }

  std::sort(event_ordering.begin(), event_ordering.end(),
            [this](const int i, const int j) {
              return this->profile[i].time < this->profile[j].time;
            });
  auto pr{List<Timepoint<T>>::tail};
  for (auto e : event_ordering) {
    profile.add_after(pr, e);
    pr = e;
  }
    for(auto i : lct_order) {
        inLeftCut[i] = true;
    }
}

template <typename T> void CumulativeEdgeFinding<T>::propagate() {

  horizontallyElasticEdgeFinderForward();

}

template <typename T>
void CumulativeEdgeFinding<T>::addTask(const int i) {
    
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "  * add " << (i) << std::endl;
  }
#endif

//    if(not checkLeftCut()) {
//        std::cout << profile;
//        assert(false);
//    }

  assert(not inLeftCut[i]);
  profile.re_add();
  profile.re_add();
  profile.re_add();

  inLeftCut[i] = true;
    
    assert(checkLeftCut());

//    if(not checkLeftCut()) {
//        std::cout << profile;
//        assert(false);
//    }
}

template <typename T> void CumulativeEdgeFinding<T>::leftCut(const int i) {

#ifdef DBG_SEF
  if (DBG_SEF) {
    if (lc_ptr != lct_order.rend())
      std::cout << "leftcut = " << *lc_ptr << ":\n" << profile;
    else
      std::cout << "leftcut = empty:\n" << profile;
  }
#endif

  if (inLeftCut[i]) {
    // task j is in the current profile, remove the tasks in (j..*ti]
    for (; *lc_ptr != i; ++lc_ptr) {
      rmTask(*lc_ptr);
    }
  } else {
    do {
      addTask(*(--lc_ptr));
    } while (*lc_ptr != i);
  }
}

template <typename T>
void CumulativeEdgeFinding<T>::rmTask(const int i) {
    
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "  * rm " << (i) << std::endl;
  }
#endif

    assert(inLeftCut[i]);

    assert(checkLeftCut());

    //    std::cout << "herea\n";
    profile.remove(lct_[i]);
    //    std::cout << "hereb\n";
    profile.remove(ect_[i]);
    //    std::cout << "herec\n";
    profile.remove(est_[i]);
    //    std::cout << "hered\n";

    inLeftCut[i] = false;

    assert(checkLeftCut());
}


template <typename T>
void CumulativeEdgeFinding<T>::addPrime(const int i) {
    
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "  * add " << i << "'\n";
  }
#endif

    auto _est{est(i)};
    auto _lct{std::min(profile.rbegin()->time, ect(i))};

    auto maxcap{capacity.max(m_solver)};
    auto ip{the_tasks.size()};

    profile[est_[ip]].capacity = profile[ect_[ip]].capacity = profile[lct_[ip]].capacity = maxcap;
    profile[est_[ip]].time = _est;
    profile[ect_[ip]].time = profile[lct_[ip]].time = _lct;
    profile[est_[ip]].increment = profile[est_[ip]].incrementMax = mindemand(i);
    profile[ect_[ip]].increment = profile[lct_[ip]].incrementMax = -mindemand(i);

    auto p{profile.rbegin()};
    while (p != profile.rend()) {
      if (p->time <= _lct) {
        profile.add_after(p.index, lct_[ip]);
        profile.add_after(p.index, ect_[ip]);
        break;
      }
      ++p;
    }
    if(p==profile.rend()) {
        profile.add_front(lct_[ip]);
        profile.add_front(ect_[ip]);
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

    assert(checkLeftCut());
}

template <typename T> void CumulativeEdgeFinding<T>::rmPrime() {

  //    std::cout << "here1\n";
  profile.remove_and_forget(est_.back());
  //    std::cout << "here2\n";
  profile.remove_and_forget(ect_.back());
  //    std::cout << "here3\n";
  profile.remove_and_forget(lct_.back());
  //    std::cout << "here4\n";

  assert(checkLeftCut());
}

template <typename T>
void CumulativeEdgeFinding<T>::horizontallyElasticEdgeFinderForward() {
    
    in_conflict.clear();

    std::sort(lct_order.begin(), lct_order.end(),
              [this](const int i, const int j) { return lct(i) < lct(j); });

    lc_ptr = lct_order.rbegin();
    
    if(lct(*(lct_order.begin())) == Constant::Infinity<T>)
        return;

  #ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "\n\nstart propagation (capacity = " << capacity.max(m_solver) << ")\n";
      for (auto j : lct_order) {
          std::cout << "task t" << std::setw(3) << std::left << the_tasks[j].id() << ": " << asciiArt(j) << std::endl;
      }
    }
  #endif

    buildFullProfile();
    
    forwardDetection();
    forwardAdjustment();
    
    
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "end propagation\n";
  }
#endif
}

template <typename T>
void CumulativeEdgeFinding<T>::forwardAdjustment() {
    
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << " adjustments on " << in_conflict << "\n";
  }
#endif

    while(not in_conflict.empty()) {
        
        auto i{in_conflict.front()};
        auto j{prec[i]};

        leftCut(j);

        scheduleOmega(capacity.max(m_solver) - mindemand(i), lct(j));

        auto t{profile.begin()->time + ceil_division(overflow, mindemand(i))};

#ifdef DBG_SEF
        if (DBG_SEF) {
          std::cout << "\n ==> adjust est(" << i << ") to " << t << std::endl;
        }
#endif

        m_solver.set(the_tasks[i].start.after(t), {this, Constant::NoHint});

        in_conflict.pop_front();
    }
}

template <typename T> void CumulativeEdgeFinding<T>::forwardDetection() {

  auto cap{capacity.max(m_solver)};
  auto stop{lct_order.rend()};
  --stop;

  // explore the tasks by decreasing lct
  for (auto ii{lct_order.rbegin()}; ii != stop;) {
    auto is{ii};
        
        
#ifdef DBG_SEF
        if (DBG_SEF) {
          std::cout << " - analyse tasks whose lct is " << lct(*ii)
                    << std::endl;
        }
#endif
        
        // remove all tasks whose lct is equal the current max
        do {

            rmTask(*ii);
            ++ii;
            ++lc_ptr;

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

              addPrime(i);

              auto omega_ect{scheduleOmega(capacity.max(m_solver), lct_j)};
#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout << " ect^H = " << omega_ect
                            << " / lct(S) = " << lct_j << std::endl;
                }
#endif

                if (omega_ect > lct_j) {

#ifdef DBG_SEF
                  if (DBG_SEF) {
                    std::cout << " task " << i
                              << " is in conflict with task interval "
                              << *(lct_order.begin()) << ".." << j << std::endl;
                  }
#endif
                  prec[i] = j;
                  in_conflict.add(i);
                } else {

#ifdef DBG_SEF
                  if (DBG_SEF) {
                    std::cout << " compute bound\n";
                  }
#endif

                  auto ti{ii};
                  computeBound(i);
                  if (beta != -1) {
#ifdef DBG_SEF
                    if (DBG_SEF) {
                      std::cout << "  - beta = " << beta << std::endl;
                    }
#endif
                    for (; *ti != beta; ++ti) {
                      rmTask(*ti);
                      ++lc_ptr;
                    }
                    auto ect_i_H{
                        scheduleOmega(capacity.max(m_solver), lct(beta))};

                    if (ect_i_H > lct(beta)) {

#ifdef DBG_SEF
                      if (DBG_SEF) {
                        std::cout << " task " << i
                                  << " is in conflict with task interval "
                                  << *(lct_order.begin()) << ".." << beta
                                  << std::endl;
                      }
#endif
                      prec[i] = beta;
                      in_conflict.add(i);
                    }
                  }
                  if (prec[i] == -1 and alpha != -1) {

#ifdef DBG_SEF
                    if (DBG_SEF) {
                      std::cout << "  - alpha = " << alpha << std::endl;
                    }
#endif
                    for (; *ti != alpha; ++ti) {
                      rmTask(*ti);
                      ++lc_ptr;
                    }
                    auto ect_i_H{
                        scheduleOmega(capacity.max(m_solver), lct(alpha))};

                    if (ect_i_H > lct(alpha)) {

#ifdef DBG_SEF
                      if (DBG_SEF) {
                        std::cout << " task " << i
                                  << " is in conflict with task interval "
                                  << *(lct_order.begin()) << ".." << alpha
                                  << std::endl;
                      }
#endif

                      prec[i] = alpha;
                      in_conflict.add(i);
                    }
                  }

                  for (; ti != ii;) {
                    --ti;
                    --lc_ptr;
                    addTask(*ti);
                  }
                }

#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout << "  * rm " << i << "'\n";
                }
#endif

                rmPrime();

                for (auto p{profile.begin()}; p != profile.end(); ++p) {
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

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmegaAdjustment(int u, const T max_lct) {
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "[schedule tasks for adjustment until " << *lc_ptr
              << "] profile=" << std::endl
              << profile;
  }
  verify();
#endif

  auto saved_size{profile.size()};
  auto sentinel{profile.end()};
  auto next{profile.begin()};

  overflow = 0;
  T omega_ect{-Constant::Infinity<T>};
  //  T S{0};
  T h_req{0};
  T hmaxInc{0};
  T overlap{0};
  T slackUnder{0};
  T slackOver{0};
  T available{0};
  T C{capacity.max(m_solver)};
  T h{mindemand(u)};

  while (next->time < max_lct) { // next != sentinel) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
    }
#endif

    t->overflow = overflow;
    auto l = next->time - t->time;

    // hmaxInc is the sum of demands of the tasks that could be processed at
    // time
    hmaxInc += t->incrementMax;
    // h_max is the min between the resource's capacity and the total of the
    // demands
    auto h_max{std::min(hmaxInc, C)};
    // h_req is the total demand counting tasks processed at their earliest
    h_req += t->increment;
    t->overlap = overlap;
    t->slackUnder = slackUnder;
    t->slackOver = slackOver;
    t->available = available;
    // h_cons is the amount of resource actually used in the optimistic scenario
    // (min between what is required + due from earlier, and what is available)
    auto h_cons{std::min(h_req + overflow, h_max)};

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << " h_max=" << h_max << ", h_req=" << h_req
                << ", h_cons=" << h_cons << ", ov=" << overflow;
      if (overflow > 0) {
        std::cout << "->" << overflow - ((h_cons - h_req) * l)
                  << " @t=" << next->time;
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
    if (h_cons > C - h) {
      t->contact = true;
    }

    // once there
    t->capacity = C - h_cons;
    overlap += std::max(h_cons - (C - h), 0) * l;
    slackOver += std::max(h_max - std::max(C - h, h_cons), 0) * l;
    slackUnder += std::max(std::min(C - h, h_max) - h_cons, 0) * l;
    available += std::min(C - h_cons, h) * l;

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

  next->overlap = overlap;
  next->slackUnder = slackUnder;
  next->slackOver = slackOver;
  next->available = available;
  if (overflow > 0)
    omega_ect = Constant::Infinity<T>;

  while (profile.size() > saved_size) {
    profile.pop_back();
  }

  return omega_ect;
}

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmega(const T C, const T max_lct) {

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "[schedule tasks until " << *lc_ptr
              << "] profile=" << std::endl
              << profile;
  }
    verify();
#endif
    
//

  auto saved_size{profile.size()};
  //  auto sentinel{profile.end()};
  //  --sentinel;

  auto next{profile.begin()};
  overflow = 0;
  T omega_ect{-Constant::Infinity<T>};
  T S{0};
  T h_req{0};

  while (next->time < max_lct) { // and next != sentinel) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
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
                << ", h_cons=" << h_cons << ", ov=" << overflow;
      if (overflow > 0) {
        std::cout << "->" << overflow - ((h_cons - h_req) * l)
                  << " @t=" << next->time;
      }
    }
#endif

    // there is some overflow, and it will be resorbed by the next time point
    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
      // then we create a new time point for the moment it will be resorbed
      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
      // h_req)));
      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
      auto new_event{//          profile.create_element(t->time + l,
                     //          t->capacity, 0, 0, 0, 0, 0, 0)};
                     profile.create_element(t->time + l, t->capacity, 0, 0)};
      profile.add_after(t.index, new_event);
      next = profile.at(new_event);
    }
    // overflow is the deficit on resource for that period (because tasks are
    // set to their earliest)profile[est_[ip]].time = _est;
    overflow += (h_req - h_cons) * l;
    // once there
    t->capacity = C - h_cons;
    //    if (overflow > 0)
    //      omega_ect = Constant::Infinity<T>;
    //    else
    if (t->capacity < C)
      omega_ect = profile[profile.next(t.index)].time;

#ifdef DBG_SEF
    if (DBG_SEF) {
      if (omega_ect != -Constant::Infinity<T>)
        std::cout << ", ect=" << omega_ect;
      std::cout << std::endl;
    }
#endif
  }

  next->overflow = overflow;
  if (overflow > 0)
    omega_ect = Constant::Infinity<T>;

  while (profile.size() > saved_size) {
    profile.pop_back();
  }

  return omega_ect;
}

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
