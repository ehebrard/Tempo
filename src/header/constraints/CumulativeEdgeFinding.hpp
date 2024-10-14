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
  Timepoint(T time, T capacity, T increment, T incrementMax, T overflow,
            T consumption, T overlap, T slackUnder, T available)
      : time(time), capacity(capacity), increment(increment),
        incrementMax(incrementMax), overflow(overflow),
        consumption(consumption), overlap(overlap), slackUnder(slackUnder),
        available(available) {}

  T time{0};
  T capacity{0};
  T increment{0};
  T incrementMax{0};
  T overflow{0};
  T consumption{0};
  T overlap{0};
  T slackUnder{0};
  T available{0};

  std::ostream &display(std::ostream &os) const {
    os << "(t=" << time << "|c=" << capacity << "|i=" << increment
       << "|im=" << incrementMax << "|of=" << overflow << "|cs=" << consumption
       << "|ol=" << overlap << "|su=" << slackUnder << "|av=" << available
       << ")";
    return os;
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Timepoint<T> &x) {
  return x.display(os);
}

template<typename T>
struct ExplanationData {
  std::vector<int> omega;
  T lb_omega;
  T ub_omega;
  int i;
  T lb_i;
};


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

  std::vector<int> prec; // mapping from the tasks that need to be adjusted to
  // the corresponding left-cut interval
  std::vector<int>
      contact; // mapping the lower bound of the set explaining the adjustment
  std::vector<T> contact_time; // mapping the lower bound of the set explaining
  // the adjustment
  SparseSet<> in_conflict; // tasks that have been found to be in conflict with
  // an left-cut interval
  std::vector<bool> inLeftCut;
  std::vector<int>::reverse_iterator
      lc_ptr; // right side of the left-cut interval
  std::vector<T> tp_attributes_est_i;
  std::vector<T> tp_attributes_ect_i;
  // helpers
  List<Timepoint<T>> profile;
  std::vector<int> est_; // index of pointer in profile
  std::vector<int> ect_;
  std::vector<int> lct_;
  // std::vector<int> lst_;

  std::vector<int> event_ordering;
  std::vector<int> lct_order;
  //  int sentinel;
  int alpha;
  int beta;

  //  std::vector<std::vector<int>> explanation_task;
  //  std::vector<T> explanation_lb;
  //  std::vector<T> explanation_ub;
  //  std::vector<int> explanation_i;
  //  std::vector<T> explanation_contact;

  std::vector<ExplanationData<T>> explanation;

  Reversible<size_t> num_explanations;

  //    int start_interval;

  // tools for adjustments
  std::vector<int>
      minEct; // the minimum ect among tasks of the set use for explanation
  std::vector<int> maxOverflow; // the number of energy units to schedule on the
                                // upper part of the profile
  std::vector<bool> isFeasible; // use to check if the correctonding scheduling
                                // is feasible (without preemption)

  T overflow;

  std::vector<std::vector<int>> triggers;

public:
  static int est_flag;
  static int ect_flag;
  static int lct_flag;
  // static int lst_flag;
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
  bool hasFixedPart(const unsigned i) const;
  T overlapedFixedPartEnergy(const unsigned i, const unsigned j) const;
  std::vector<T> overlapedFixedPartEnergy();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;
  //  T scheduleOmega(const T C, const T max_lct);
  void computeBound(const int i);

  /* Tools for the adjustments */
  T scheduleOmegaDetection(
      const T C, const int i
//      ,const T max_lct
                           ); // scheduling omega used for the detection
  T scheduleOmegaAdjustment(
      const T C, const int i
//      ,const T max_lct
                            ); // scheduling omega used for the adjustment
  int computeEstPrime(
      const int i); // end of the scheduling in the upper part of the profile
  void computeMaximumOverflow(
      const int i); // energy to spend at the upper part of the profile
  std::vector<T>
  consOverSlackAvail(const int time); // retrive the attributes of the time
                                      // points not belong to the profile

  void buildFullProfile();

  void addTask(const int i);
  void rmTask(const int i);

  void leftCut(const int i);

  bool addPrime(const int i);
  void rmPrime();



  bool addTimePoint_i(const int i);
  void rmTimePoint_i();

  void horizontallyElasticEdgeFinderForward();
  void forwardDetection();
  void forwardAdjustment();
    
    void overloadBound();

  // function used in explanation
  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;
  std::string asciiArt(const int i) const;

  std::vector<bool> inprof;
  bool verify() {
    return true;
    inprof.clear();
    inprof.resize(the_tasks.size() + 4, false);
    auto prev_t{-Constant::Infinity<T>};
    for (auto p{profile.begin()}; p != profile.end(); ++p) {
      if (inprof[p.index]) {
        std::cout << p.index << " is in the list twice (loop)\n";
        std::cout << profile << std::endl;
        return false;
      }

      inprof[p.index] = true;

      if (prev_t > p->time) {
        std::cout << p.index << " is wrongly ordered\n";
        return false;
      }
    }

    return checkLeftCut();
  }

  bool checkLeftCut() {
    int count{0};
    for (auto p{profile.begin()}; p != profile.end(); ++p) {
      if (/*p.index != sentinel and*/ p.index != est_.back() and
          p.index != ect_.back() and p.index != lct_.back())
        ++count;

      //        std::cout << "count " << p.index << " (" << *p << " / " <<
      //        sentinel << ")\n";
    }
    bool ok{true};
    for (unsigned i{0}; ok and i < the_tasks.size(); ++i) {
      if (inLeftCut[i]) {
        count -= 3;
        if (not(inprof[est_[i]] and inprof[ect_[i]] and inprof[lct_[i]])) {

          std::cout << "not all ptrs of task " << the_tasks[i]
                    << " are in the profile (est[" << i
                    << "]:" << inprof[est_[i]] << ", ect[" << i
                    << "]:" << inprof[ect_[i]] << ", lct[" << i
                    << "]:" << inprof[lct_[i]] << ")\n";

          ok = false;
        }
      }
    }

    if (count != 0) {

      std::cout << "there are " << count
                << " too many elements in the profile\n";

      return false;
    }

    //    ok &= (count == 0);
    return ok;
  }
};






template <typename T>
int CumulativeEdgeFinding<T>::est_flag = 0;

template <typename T>
int CumulativeEdgeFinding<T>::ect_flag = 1;

template <typename T>
int CumulativeEdgeFinding<T>::lct_flag = 2;

/*template <typename T>
int CumulativeEdgeFinding<T>::lst_flag = 3;*/

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
  if (lct(i) == Constant::Infinity<T>) {
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
bool CumulativeEdgeFinding<T>::hasFixedPart(const unsigned i) const {
  return lst(i) < ect(i);
}

template <typename T>
T CumulativeEdgeFinding<T>::overlapedFixedPartEnergy(const unsigned i, const unsigned j) const {
  auto energy{0};
  if (lct(i) > lct(j) and hasFixedPart(i) and lst(i) < lct(j)) {
    energy = (std::max(0, std::min(ect(i), lct(j)) - lst(i))) * mindemand(i);
  }
  return energy;
}

template <typename T>
std::vector<T> CumulativeEdgeFinding<T>::overlapedFixedPartEnergy() {
  std::vector<T> energy;
  energy.resize(the_tasks.size());
  for (auto j : lct_order) {
    auto E{0};
    for (auto i = j + 1; i < the_tasks.size(); ++i) {
      E += overlapedFixedPartEnergy(i, j);
    }
    energy[j] = E;
    while (j + 1 < the_tasks.size() and lct(j + 1) == lct(j)) {
      energy[j + 1] = energy[j];
      ++j;
    }
  }
  return energy;
}

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeEdgeFinding<T>::CumulativeEdgeFinding(
    Solver<T> &solver, const Interval<T> sched, const NumericVar<T> cap,
    const ItTask beg_task, const ItTask end_task, const ItNVar beg_dem,
    const ItBVar beg_disj)
    : m_solver(solver), num_explanations(0, &(m_solver.getEnv())) {
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
  // lst_.resize(ip + 1);

  contact.resize(ip);
  contact_time.resize(ip);
  minEct.resize(ip);
  isFeasible.resize(ip);
  maxOverflow.resize(ip);
  tp_attributes_est_i.resize(3);
  tp_attributes_ect_i.resize(3);

  auto maxcap{capacity.max(m_solver)};
  for (unsigned i = 0; i < ip; ++i) {
    precedence[i].resize(the_tasks.size());
  }
  for (unsigned i = 0; i < ip; ++i) {
    precedence[i].resize(the_tasks.size());
    est_[i] = profile.create_element(est(i), maxcap, mindemand(i), mindemand(i),
                                     0, 0, 0, 0, 0);
    ect_[i] =
        profile.create_element(ect(i), maxcap, -mindemand(i), 0, 0, 0, 0, 0, 0);
    lct_[i] =
        profile.create_element(lct(i), maxcap, 0, -mindemand(i), 0, 0, 0, 0, 0);
    /*if (hasFixedPart(i))
    lst_[i] =
    profile.create_element(lst(i), maxcap, mindemand(i), mindemand(i), 0, 0, 0,
    0, 0);*/
  }
  //  sentinel = profile.create_element(Constant::Infinity<T>, 0, 0, 0, 0, 0, 0,
  //  0, 0);

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
  est_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0, 0);
  ect_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0, 0);
  lct_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0, 0);
  // lst_[ip] = profile.create_element(0, 0, 0, 0, 0, 0, 0, 0, 0);
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
  /*for (size_t i{0}; i < the_tasks.size(); ++i) {
  auto k{m_solver.wake_me_on(ub<T>(the_tasks[i].start.id()), this->id())};
  if (triggers.size() <= k)
  triggers.resize(k + 1);
  triggers[k].push_back(lst_flag + 4 * i);
  }*/
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
    /*if (hasFixedPart(i))
    profile[lst_[i]].time = lst(i);*/
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
  for (auto i : lct_order) {
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
  assert(not inLeftCut[i]);
#endif

  profile.re_add();
  profile.re_add();
  profile.re_add();

  inLeftCut[i] = true;

#ifdef DBG_SEF
  assert(verify());
#endif
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

  assert(inLeftCut[i]);
  if (not verify()) {
    std::cout << profile << std::endl;
    for (unsigned k{0}; k < inLeftCut.size(); ++k) {
      std::cout << the_tasks[k] << " " << est_[k] << "/" << ect_[k] << "/"
                << lct_[k] << std::endl;
    }
    exit(1);
  }
#endif

  //    std::cout << "herea\n";
  profile.remove(lct_[i]);
  //    std::cout << "hereb\n";
  profile.remove(ect_[i]);
  //    std::cout << "herec\n";
  profile.remove(est_[i]);
  //    std::cout << "hered\n";

  inLeftCut[i] = false;

#ifdef DBG_SEF
  assert(verify());
#endif
}

template <typename T> bool CumulativeEdgeFinding<T>::addTimePoint_i(const int i) {

  auto _est{est(i)};
  auto _ect{std::min(lct(prec[i]), ect(i))};

  if (_est >= lct(prec[i])) {

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "  * ignore " << i << "''\n";
    }
#endif

    return false;
  }

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "  * add " << i << "''\n";
  }
#endif

  auto maxcap{capacity.max(m_solver)};
  auto ip{the_tasks.size()};

  profile[est_[ip]].capacity = profile[ect_[ip]].capacity = maxcap;
  profile[est_[ip]].time = _est;
  profile[ect_[ip]].time = _ect;
  profile[est_[ip]].increment = profile[est_[ip]].incrementMax = 0;
  profile[ect_[ip]].increment = profile[ect_[ip]].incrementMax = 0;

  auto p{profile.rbegin()};
  while (p != profile.rend()) {
    if (p->time <= _ect) {
      profile.add_after(p.index, ect_[ip]);
      break;
    }
    ++p;
  }
  if (p == profile.rend()) {
    profile.add_front(ect_[ip]);
  }
  while (p != profile.rend()) {
    if (p->time <= _est) {
      profile.add_after(p.index, est_[ip]);
      break;
    }
    ++p;
  }
  if (p == profile.rend()) {
    profile.add_front(est_[ip]);
  }

#ifdef DBG_SEF
  assert(verify());
#endif

  return true;
}



template <typename T> bool CumulativeEdgeFinding<T>::addPrime(const int i) {

  auto _est{est(i)};
  auto _lct{std::min(profile.rbegin()->time, ect(i))};

  if (_est >= _lct) {

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "  * ignore " << i << "'\n";
    }
#endif

    return false;
  }

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "  * add " << i << "'\n";
  }
#endif

  auto maxcap{capacity.max(m_solver)};
  auto ip{the_tasks.size()};

  profile[est_[ip]].capacity = profile[ect_[ip]].capacity =
      profile[lct_[ip]].capacity = maxcap;
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
  if (p == profile.rend()) {
    profile.add_front(lct_[ip]);
    profile.add_front(ect_[ip]);
  }
  while (p != profile.rend()) {
    if (p->time <= _est) {
      profile.add_after(p.index, est_[ip]);
      break;
    }
    ++p;
  }
  if (p == profile.rend()) {
    profile.add_front(est_[ip]);
  }

#ifdef DBG_SEF
  assert(verify());
#endif

  return true;
}

template <typename T> void CumulativeEdgeFinding<T>::rmPrime() {

  //    std::cout << "here1\n";
  profile.remove_and_forget(est_.back());
  //    std::cout << "here2\n";
  profile.remove_and_forget(ect_.back());
  //    std::cout << "here3\n";
  profile.remove_and_forget(lct_.back());
  //    std::cout << "here4\n";

#ifdef DBG_SEF
  assert(verify());
#endif
}

template <typename T> void CumulativeEdgeFinding<T>::rmTimePoint_i() {

  //    std::cout << "here1\n";
  profile.remove_and_forget(est_.back());
  //    std::cout << "here2\n";
  profile.remove_and_forget(ect_.back());
  //    std::cout << "here3\n";
  //profile.remove_and_forget(lct_.back());
  //    std::cout << "here4\n";

#ifdef DBG_SEF
  assert(verify());
#endif
}

template <typename T>
void CumulativeEdgeFinding<T>::horizontallyElasticEdgeFinderForward() {

  in_conflict.clear();

  std::sort(lct_order.begin(), lct_order.end(),
            [this](const int i, const int j) { return lct(i) < lct(j); });

  lc_ptr = lct_order.rbegin();

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "\n\nstart propagation\n";
    for (auto j : lct_order) {
      std::cout << "task " << j << ": " << asciiArt(j) << std::endl;
    }
  }
#endif

  buildFullProfile();
    
    overloadBound();
    
  //std::cout<< " detection  " << std::endl;
  forwardDetection();
  //std::cout<< " adjustment  " << std::endl;
  forwardAdjustment();
}

template <typename T>
std::vector<T> CumulativeEdgeFinding<T>::consOverSlackAvail(const int time) {
  std::vector<T> attributes;
  attributes.reserve(4);
  auto next{profile.begin()};
  while (next != profile.end()) {
    auto t{next};
    ++next;
    if (t->time <= time && next->time > time) {
      break;
    }
  }

  --next;
    
    //std::cout << "consOverSlackAvail: " << *next << std::endl;
    
    
  attributes.push_back(next->consumption);
  attributes.push_back(next->overlap);
  attributes.push_back(next->slackUnder);
  attributes.push_back(next->available);
  return attributes;
}

template <typename T>
void CumulativeEdgeFinding<T>::forwardAdjustment() {

  // lb_expl.clear();
  // ub_expl.clear();
  // pruned_task_expl.clear();
//  auto ip{static_cast<int>(the_tasks.size())};

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << " adjustments on " << in_conflict << "\n";
  }
#endif
  while (not in_conflict.empty()) {
    auto i{in_conflict.front()};
      
#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << " adjustment on " << i << "\n";
  }
#endif
      
    auto j{prec[i]};
    auto estPrime{est(i)};
    bool pruning_flag{false};
    leftCut(j);
    if (addTimePoint_i(i))
    {
//        auto ect_h{
            scheduleOmegaAdjustment(capacity.max(m_solver), i
                                    //, lct(j)
                                    )
//        }
        ;
      computeMaximumOverflow(i);
      /*if (profile[est_[ip]].time != est(i) )
        std::cout<< profile[est_[ip]] << " contradiction  ++++++++++++++++++++ " << est(i) << " contact  " << contact_time[i] << std::endl;
      if (profile[ect_[ip]].time <= lct(j))
        std::cout<< profile[ect_[ip]] << " contradiction  ++++++++++++++++++++ " << ect(i) << " lct_j " << profile[lct_[j]] << std::endl;*/
      estPrime = computeEstPrime(i);
      if (estPrime > est(i)) {
        pruning_flag = true;

        m_solver.set(the_tasks[i].start.after(est(i)+1),
                     {this, static_cast<hint>(explanation.size() - 1)});
      }
      rmTimePoint_i();
    }


    /*leftCut(j);
    bool pruning_flag{false};
    if (addTimePoint_i(i))
    {
      auto ect_h{scheduleOmegaAdjustment(capacity.max(m_solver), i, lct(j))};
      computeMaximumOverflow(i);
      auto estPrime{est(i)};
      if (maxOverflow[i] > 0) {
        estPrime = computeEstPrime(i);
        if (estPrime > est(i)) {

#ifdef DBG_SEF
          if (DBG_SEF) {
            std::cout << "\n ==> adjust est(" << the_tasks[i].id() << ") to "
                      << estPrime << " (was " << est(i) << ") :: "
                      << m_solver.pretty(the_tasks[i].start.after(estPrime))
                      << std::endl;
          }
#endif

          pruning_flag = true;

          m_solver.set(the_tasks[i].start.after(estPrime),
                       {this, static_cast<hint>(explanation.size() - 1)});
        }
      }
      rmTimePoint_i();
    }*/
    if (not pruning_flag) {
      explanation.pop_back();
      //      explanation_task.pop_back();
      //      explanation_lb.pop_back();
      //      explanation_ub.pop_back();
      //        explanation_i.pop_back();
      //        explanation_contact.pop_back();
    } else {
      ++num_explanations;
    }

    in_conflict.pop_front();
  }
}

template <typename T>
int CumulativeEdgeFinding<T>::computeEstPrime(const int i) {
  //std::cout<< " compute est prime " << std::endl;
  auto cont{contact_time[i]};
  auto Ov{maxOverflow[i]};
  auto ip{static_cast<int>(the_tasks.size())};
  auto next{profile.at(contact[i])};
  if (cont < est(i))
    next = profile.at(est_[ip]);
  auto estPrime{-Constant::Infinity<T>};
  while (next != profile.end()) {
    auto t{next};
    ++next;
    auto overl{next->overlap - t->overlap};
    assert(overl >= 0);
    //std::cout<< " ------> compute est prime: overlap " << overl  << " Ov " << Ov << " consomption " << (next->consumption)
    //<< " time point " << *next << " C-h " <<  ((capacity.max(m_solver) - mindemand(i)))<<  std::endl;
    if (Ov > overl) {
      assert(Ov >= overl);
      Ov -= overl;
    } else {
      /*std::cout<< " compute est prime " << (ceil_division(Ov, (t->consumption -
                                                       (capacity.max(m_solver) -
                                                        mindemand(i))))) << " consomption " << (t->consumption) << std::endl;*/
      /*if (t->consumption <= capacity.max(m_solver) - mindemand(i))
      {
        std::cout<< "Task_i : est = " << est(i) << " -- lct = " << lct(i) << " -- d_i = " << minduration(i)
        << " -- c_i = " << mindemand(i) << " -- " << " Capa " << capacity.max(m_solver) <<
          " max Ov " << maxOverflow[i] <<   std::endl;
        for (auto j:lct_order)
        {
          if (ect(j) >= contact_time[i] and lct(j) <= lct(prec[i]))
          {
            std::cout<< " est = " << est(j) << " -- lct = " << lct(j) << " -- d_i = " << minduration(j) << " -- c_i = " << mindemand(j) << std::endl;
          }
        }
        for (auto tp:profile)
        {
          if (tp.time >= 20 and tp.time <= 34)
            std::cout<< tp << std::endl;
        }


      }*/
      if (t->consumption > capacity.max(m_solver) -mindemand(i))
      {
        estPrime = std::min(next->time,
                           t->time + ceil_division(Ov, (t->consumption -
                                                        (capacity.max(m_solver) -
                                                         mindemand(i)))));
        return estPrime;
      }
    }
  }
  return estPrime;
}

template <typename T> void CumulativeEdgeFinding<T>::overloadBound() {
    
#ifdef DBG_SEF
    if (DBG_SEF) {
    std::cout << "compute overload bound\n";
        for (auto j : lct_order) {
          std::cout << "task " << j << ": " << asciiArt(j) << std::endl;
        }
    }
#endif
    
    auto cap{capacity.max(m_solver)};
    auto last{*lct_order.rend()};
    auto omega_ect = scheduleOmegaDetection(cap, last
//                                            , lct(last)
                                            );
    
    if(omega_ect == Constant::Infinity<T>) {
        
        
        std::cout << "overload fail!\n";
        
        throw Failure<T>({this, Constant::NoHint});
    } else if(omega_ect > schedule.end.max(m_solver)) {
        
        std::cout << "makespan pruning!\n";
        
        m_solver.set(schedule.end.after(omega_ect), {this, Constant::NoHint});
    } else {
#ifdef DBG_SEF
        if (DBG_SEF) {
            std::cout << "nothing\n";
        }
#endif
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
      std::cout << " - analyse tasks whose lct is " << lct(*ii) << std::endl;
    }
#endif

    // remove all tasks whose lct is equal the current max
    do {

      rmTask(*ii);
      ++ii;
      ++lc_ptr;

    } while (ii != lct_order.rend() and lct(*is) == lct(*ii));

    // if there are no more tasks, all those in the current level have the same
    // lct and we can stop
    if (ii == lct_order.rend())
      break;

    // lct of the task that precedes all the task whose lct is lct(i)
    auto j{*ii};
    auto lct_j{lct(j)};

    // otherwise, add their "prime" versions one by one and run scheduleOmega
    while (is != ii) {

      auto i{*is};

      if (ect(i) < lct(i)) {
        // std::cout<< profile[est_[the_tasks.size()]] << std::endl;
        for (auto p{profile.begin()}; p != profile.end(); ++p) {
          p->capacity = cap;
        }
        auto omega_ect{0};
        if (addPrime(i))
        {
          omega_ect = scheduleOmegaDetection(cap, i
//                                             , lct_j
                                             );
          rmPrime();
        }
          /*if (omega_ect > lct(j))
          {
            std::cout<< profile[est_[the_tasks.size()]] << " --- " << profile[ect_[the_tasks.size()]] << " --- " << profile[lct_[the_tasks.size()]] << std::endl;
            std::cout<< est(i) << " --- " << lct(i) << " --- " << minduration(i) << std::endl;
          }*/


#ifdef DBG_SEF
          if (DBG_SEF) {
            std::cout << " ect^H = " << omega_ect << " / lct(S) = " << lct_j
                      << std::endl;
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
            //rmPrime();
          }

          else {

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
              for (auto p{profile.begin()}; p != profile.end(); ++p) {
                p->capacity = cap;
              }
              auto ect_i_H{0};
              if (addPrime(i))
              {
                ect_i_H = scheduleOmegaDetection(cap, i
                                                 //, lct(beta)
                                                 );
                rmPrime();
              }

              if (ect_i_H > lct(beta)) {
#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout
                      << " task " << i << " is in conflict with task interval "
                      << *(lct_order.begin()) << ".." << beta << std::endl;
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
              for (auto p{profile.begin()}; p != profile.end(); ++p) {
                p->capacity = cap;
              }
              auto ect_i_H{0};
              if (addPrime(i))
              {
                ect_i_H = scheduleOmegaDetection(cap, i
//                                                 , lct(alpha)
                                                 );
                rmPrime();
              }


              if (ect_i_H > lct(alpha)) {

#ifdef DBG_SEF
                if (DBG_SEF) {
                  std::cout
                      << " task " << i << " is in conflict with task interval "
                      << *(lct_order.begin()) << ".." << alpha << std::endl;
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

          /*for (auto p{profile.begin()}; p != profile.end(); ++p) {
            p->capacity = cap;
          }*/
        //}
      }

      ++is;
    }
  }
}

template <typename T>
void CumulativeEdgeFinding<T>::computeMaximumOverflow(const int i) {
  auto e_idx{num_explanations};
  if (explanation.size() <= e_idx) {
    explanation.resize(e_idx + 1);
  } else {
    explanation[e_idx].omega.clear();
  }
  auto ECT{Constant::Infinity<T>};
  auto minEnergy{Constant::Infinity<T>};
  auto lb{Constant::Infinity<T>};
  auto k{prec[i]};

  for (auto j : lct_order) {
    if (lct(j) > lct(k))
      break;
    if (ect(j) >= contact_time[i] and lct(j) <= lct(k)) {
      ECT = std::min(ect(j), ECT);
      minEnergy = std::min(energy(j), minEnergy);
      lb = std::min(lb, est(j));
      explanation[e_idx].omega.push_back(j);
    }
  }

  explanation[e_idx].lb_omega = lb;
  explanation[e_idx].ub_omega = lct(k);
  explanation[e_idx].i = i;
  explanation[e_idx].lb_i = std::min(est(i), contact_time[i]);

  minEct[i] = ECT;


  auto cont{contact_time[i]};
  auto ip{static_cast<int> (the_tasks.size())};

  if (ect(i) < lct(k)) {
    //std::cout << "i:" << i << " // " << cont << " <> " << est(i) << std::endl;

    /*if (cont < est(i)) {
      cont = est(i);
    }*/
    //      std::cout << "consOverSlackAvail(" << cont << ")\n";
    auto overlap1{0};
    auto slackUnder1{0};
    auto available1{0};

    if (cont > est(i))
    {
      auto t{profile[contact[i]]};
      overlap1 = t.overlap;
      slackUnder1 = t.slackUnder;
      available1 = t.available;
    } else
    {
      auto t{profile[est_[ip]]};
      overlap1 = t.overlap;
      slackUnder1 = t.slackUnder;
      available1 = t.available;
    }

    //      std::cout << "consOverSlackAvail(" << ect(i) << ")\n";
    //auto attributes2{consOverSlackAvail(ect(i))};
    //auto overlap2{tp_attributes_ect_i[0]};
    auto t{profile[ect_[ip]]};
    auto overlap2 = t.overlap;

    //      std::cout << "consOverSlackAvail(" << lct(k) << ")\n";
    //auto attributes3{consOverSlackAvail(lct(k))};
    t = profile[lct_[k]];
    auto slackUnder2{t.slackUnder};
    auto overlap3{t.overlap};
    auto available2{t.available};
    if (available2 - available1 < energy(i)) {
      if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
        maxOverflow[i] = overlap3 - overlap1;
      } else {
        maxOverflow[i] = overlap3 - overlap1 - (slackUnder2 - slackUnder1);
      }
    } else {
      if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
        maxOverflow[i] = overlap2 - overlap1;
      } else {
        maxOverflow[i] = overlap2 - overlap1 - (slackUnder2 - slackUnder1);
      }
    }
  } else {
    auto overlap1{0};
    auto slackUnder1{0};
    if (cont > est(i))
    {
      auto t{profile[contact[i]]};
      overlap1 = t.overlap;
      slackUnder1 = t.slackUnder;
    } else
    {
      auto t{profile[est_[ip]]};
      overlap1 = t.overlap;
      slackUnder1 = t.slackUnder;
    }
    auto t = profile[lct_[k]];
    auto slackUnder2{t.slackUnder};
    auto overlap2{t.overlap};
    if (isFeasible[i] and minEnergy > slackUnder2 - slackUnder1) {
      maxOverflow[i] = overlap2 - overlap1;
    } else {
      maxOverflow[i] = overlap2 - overlap1 - (slackUnder2 - slackUnder1);
    }
  }

  if (maxOverflow[i] < 0 /*or m_solver.num_cons_propagations == 70*/) {
    std::cout << "bug (" << m_solver.num_cons_propagations <<  " est_i " << est(i) <<  " ect_i " << ect(i)  << ")\n";

    std::cout << "max_Overflow_is_equal = " << maxOverflow[i]  << " contact time " << contact_time[i] << profile[contact[i]] << " --- " << profile[est_[ip]] << ")\n";

    for (auto j : lct_order) {
      if (lct(j) > lct(k))
        break;
      //if (ect(i) >= contact_time[i])
        std::cout << std::setw(3) << j << ": " << asciiArt(j) << std::endl;
    }
    std::cout << std::endl
              << std::setw(3) << i << ": " << asciiArt(i) << std::endl;

    std::cout << "capacity = " << capacity.max(m_solver)
              << " id = " << this->id() << std::endl;

    std::cout << profile << std::endl;

    exit(1);
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
    E += energy(j);
    if (lct(j) <= ect(i) and est(i) < lct(j)) {
      auto slack{(capacity.max(m_solver) - mindemand(i)) * lct(j) - E};
      if (slack < minSlack[0] and profile[lct_[j]].overflow > 0) {
        minSlack[0] = slack;
        alpha = j;
      }
    }
    // else
    if (lct(j) > ect(i)) {
      auto slack{capacity.max(m_solver) * lct(j) - E};
      if (slack < minSlack[1] and profile[lct_[j]].overflow > 0) {
        minSlack[1] = slack;
        beta = j;
      }
    }
  }
}

// template <typename T> T CumulativeEdgeFinding<T>::scheduleOmega(const T C,
// const T max_lct) {
//
//
//#ifdef DBG_SEF
//   if (DBG_SEF) {
//     std::cout << "[schedule tasks until " << *lc_ptr
//               << "] profile=" << std::endl
//               << profile;
//   }
//   assert(verify());
//#endif
//
////
//
//  auto saved_size{profile.size()};
////  auto sentinel{profile.end()};
////  --sentinel;
//
//  auto next{profile.begin()};
//  overflow = 0;
//  T omega_ect{-Constant::Infinity<T>};
//  T S{0};
//  T h_req{0};
//
//  while (next->time < max_lct) {
//    auto t{next};
//    ++next;
//
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
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
//      std::cout << " h_max=" << h_max << ", h_req=" << h_req
//                << ", h_cons=" << h_cons << ", ov=" << overflow;
//      if (overflow > 0) {
//        std::cout << "->" << overflow - ((h_cons - h_req) * l)
//                  << " @t=" << next->time;
//      }
//    }
//#endif
//
//    // there is some overflow, and it will be resorbed by the next time point
//    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
//      // then we create a new time point for the moment it will be resorbed
//      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
//      // h_req)));
//      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
//      auto new_event{
//          profile.create_element(t->time + l, t->capacity, 0, 0, 0, 0, 0, 0,
//          0)};
//      profile.add_after(t.index, new_event);
//      next = profile.at(new_event);
//    }
//    // overflow is the deficit on resource for that period (because tasks are
//    // set to their earliest)profile[est_[ip]].time = _est;
//    overflow += (h_req - h_cons) * l;
//    // once there
//    t->capacity = C - h_cons;
//    if (overflow > 0)
//      omega_ect = Constant::Infinity<T>;
//    else if (t->capacity < C)
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
//  }
//
//  return omega_ect;
//}
//
//
//

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmegaDetection(const T C, const int i
//                                                   ,const T max_lct
                                                   ) {

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "[schedule tasks until " << *lc_ptr
              << "] profile=" << std::endl
              << profile;
  }
  assert(verify());
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

  // while (next->time <= max_lct and next != ) {
  while (next != profile.end()) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
    }
#endif

    t->overflow = overflow;
    auto l = next->time - t->time;
//      if(next->time == Constant::Infinity<T>) {
//
//      }
      

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

      
      
      
//      auto absorbtion{(h_cons > h_req) and (l > 0)};
      
//      if(absorbtion)
      auto al{std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req))};
      
//      std::cout << ", l=" << l << ", al="  ;
      
    // there is some overflow, and it will be resorbed by the next time point
    if (overflow > 0 and
        al < l) {
//        overflow < ((h_cons - h_req) * l)) {
      // then we create a new time point for the moment it will be resorbed
      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
      // h_req)));
        l = al;
        
        
//#ifdef DBG_SEF
//    if (DBG_SEF) {
//      std::cout << " overflow absorbtion event @t" << t->time + l << std::endl;
//    }
//#endif
        
      auto new_event{profile.create_element(t->time + l, t->capacity, 0, 0, 0,
                                            0, 0, 0, 0)};
      profile.add_after(t.index, new_event);
      next = profile.at(new_event);
    }
      
//      std::cout << l;
      
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
      std::cout << " h_max=" << h_max << ", h_req=" << h_req
                << ", h_cons=" << h_cons
                << ", ov=" << (overflow - (h_req - h_cons) * l);
      if (overflow > 0) {
        std::cout << "->" << overflow << " @t=" << next->time;
      }

      if (omega_ect != -Constant::Infinity<T>)
        std::cout << ", ect=" << omega_ect;
      std::cout << std::endl;
    }
#endif
  }

  if (overflow > 0) {
    //      contact[i] = -1;
    //      auto t{profile.reverse_at(next.index)};
    //      while(t != )
    //
    //
    //      auto t{next};
    //      --t;
    //      while(t )

    contact[i] = -1;
    auto previous{next};
    --previous;
    while (previous != profile.begin()) {
      auto t{previous};
      --previous;
      if (t->overflow < overflow) {

        contact_time[i] = t->time;

        if (t.index > 3 * static_cast<int>(the_tasks.size()))
          --t;

        contact[i] = t.index;

        break;
      }
    }
    if (contact[i] == -1 and previous == profile.begin()) {
      contact[i] = previous.index;
      contact_time[i] = previous->time;
    }
      
      //std::cout << "contact = " << profile[contact[i]] << " contact time = " << contact_time[i] << std::endl;
  }

  while (profile.size() > saved_size) {
    profile.pop_back();
  }

  return omega_ect;
}

template <typename T>
T CumulativeEdgeFinding<T>::scheduleOmegaAdjustment(const T C, const int i
//                                                    ,const T max_lct
                                                    ) {

#ifdef DBG_SEF
  if (DBG_SEF) {
    std::cout << "[schedule tasks until " << *lc_ptr
              << "] profile=" << std::endl
              << profile;
  }
  assert(verify());
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
  T overlap{0};
  T slackUnder{0};
  T available{0};
  T h{mindemand(i)};
  isFeasible[i] = true;
  tp_attributes_est_i.resize(3);
  tp_attributes_ect_i.resize(3);

    T prev_cons{0};
    
  // while (next->time <= max_lct) {
  while (next != profile.end()) {
    auto t{next};
    ++next;

#ifdef DBG_SEF
    if (DBG_SEF) {
      std::cout << "jump to e_" << t.index << " = t_" << t->time << ":";
    }
#endif

    t->overflow = overflow;
    auto l = next->time - t->time;
    next->overlap = overlap;
    next->slackUnder = slackUnder;
    next->available = available;

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
    t->consumption = h_cons;

    // there is some overflow, and it will be resorbed by the next time point
    if (overflow > 0 and overflow < ((h_cons - h_req) * l)) {
      // then we create a new time point for the moment it will be resorbed
      // l = std::max(Gap<T>::epsilon(), ceil_division<T>(overflow, (h_cons -
      // h_req)));
      l = std::max(Gap<T>::epsilon(), overflow / (h_cons - h_req));
      auto new_event{profile.create_element(t->time + l, t->capacity, 0, 0, 0,
                                            0, 0, 0, 0)};
      profile.add_after(t.index, new_event);
      next = profile.at(new_event);
    }
    // overflow is the deficit on resource for that period (because tasks are
    // set to their earliest)profile[est_[ip]].time = _est;
    overflow += (h_req - h_cons) * l;
    // once there
    t->capacity = C - h_cons;
    // check the feasibility of the scheduling
    if (overflow > 0)
      isFeasible[i] = false;

//      if(l>0) {
          // update the values of overlap, slack under and available
      overlap += (std::max(h_cons - (C - h), 0) * l);
      slackUnder += (std::max(std::min(C - h, h_max) - h_cons, 0) * l);
      available += (std::min(C - h_cons, h) * l);
    if ((t->time) < est(i))
    {
      if (next->time == est(i))
      {
        tp_attributes_est_i[0] = overlap;
        tp_attributes_est_i[1] = slackUnder;
        tp_attributes_est_i[2] = available;
      } else if (next->time > est(i))
      {
        tp_attributes_est_i[0] = t->overlap + (std::max(h_cons - (C - h), 0) * (est(i) - t->time));
        tp_attributes_est_i[1] = t->slackUnder + (std::max(std::min(C - h, h_max) - h_cons, 0) * (est(i) - t->time));
        tp_attributes_est_i[2] = t->available + (std::min(C - h_cons, h) * (est(i) - t->time));
      }
    }

    if ((t->time) < ect(i))
    {
      if (next->time == ect(i))
      {
        tp_attributes_ect_i[0] = overlap;
        tp_attributes_ect_i[1] = slackUnder;
        tp_attributes_ect_i[2] = available;
      } else if (next->time > ect(i))
      {
        tp_attributes_ect_i[0] = t->overlap + (std::max(h_cons - (C - h), 0) * (ect(i) - t->time));
        tp_attributes_ect_i[1] = t->slackUnder + (std::max(std::min(C - h, h_max) - h_cons, 0) * (ect(i) - t->time));
        tp_attributes_ect_i[2] = t->available + (std::min(C - h_cons, h) * (ect(i) - t->time));
      }
    }
      
   
//      }
      

    if (overflow > 0)
      omega_ect = Constant::Infinity<T>;
    else if (t->capacity < C)
      omega_ect = profile[profile.next(t.index)].time;

#ifdef DBG_SEF
    if (DBG_SEF) {

      std::cout << "l=" << l << " h_max=" << h_max << ", h_req=" << h_req
                << ", h_cons=" << h_cons << "|" << prev_cons
                << ", ov=" << (overflow - (h_req - h_cons) * l);
      if (overflow > 0) {
        std::cout << "->" << overflow << " @t=" << next->time;
      }
      std::cout << ", overlap=" << overlap
                << ", slackUnder=" << slackUnder
                << ", available=" << available;

      if (omega_ect != -Constant::Infinity<T>)
        std::cout << ", ect=" << omega_ect;
      std::cout << std::endl;
    }
#endif
      
      
      if(l>0) {
          prev_cons = h_cons;
      }
          
  }
  // std::cout << " is Feasible " << " is Feasible "<< isFeasible[i] << "\n";

  while (profile.size() > saved_size) {
    profile.pop_back();
  }

  return omega_ect;
}
template <typename T>
void CumulativeEdgeFinding<T>::xplain(const Literal<T> l, const hint h,
                                      std::vector<Literal<T>> &Cl) {

  //      assert(l != Solver<T>::Contradiction);

  //  std::cout << h << "/" <<

#ifdef DBG_EXPLCE
  std::cout << "explain (" << explanation[h].size() << ") "
            << m_solver.pretty(l) << ":\n";
#endif
    
    if(l == Solver<T>::Contradiction) {
        std::cout << "xplain contradiction: TODO\n";
        exit(1);
    } else if(l.variable() == schedule.end.id()) {
        std::cout << "xplain global bound: TODO\n";
        exit(1);
    } else {
        
        for (auto i : explanation[h].omega) {
            
#ifdef DBG_EXPLCE
            std::cout << " * "
            << m_solver.pretty(the_tasks[i].start.after(explanation_lb[h]))
            << " and "
            << m_solver.pretty(the_tasks[i].end.before(explanation_ub[h]))
            << std::endl;
#endif
            
            Cl.push_back(the_tasks[i].start.after(explanation[h].lb_omega));
            Cl.push_back(the_tasks[i].end.before(explanation[h].ub_omega));
        }
        
#ifdef DBG_EXPLCE
        std::cout << " AND "
        << m_solver.pretty(
                           the_tasks[explanation[h].i].end.after(explanation[h].lb_i))
        << std::endl;
#endif
        
        Cl.push_back(the_tasks[explanation[h].i].start.after(explanation[h].lb_i));
    }

  //      Cl.push_back(geq<T>(l.variable(), explanation_lb[h]));
}

template <typename T>
std::ostream &CumulativeEdgeFinding<T>::display(std::ostream &os) const {
  os << "Cumulative Edge-Finding"; // "data/sample/j309_5.sm"

#ifdef DBG_CEDGEFINDING
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
std::ostream &CumulativeEdgeFinding<T>::print_reason(std::ostream &os, const hint) const {
  os << "cumulative-edge-finding";
  return os;
}

} // namespace tempo

#endif
