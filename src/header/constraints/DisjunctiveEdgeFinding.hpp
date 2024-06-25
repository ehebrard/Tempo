#ifndef TEMPO_DISJUNCTIVEEDGEFINDING_HPP
#define TEMPO_DISJUNCTIVEEDGEFINDING_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Solver.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"
#include "util/ThetaTree.hpp"


namespace tempo {

template <typename T> class DisjunctiveEdgeFinding : public Constraint<T> {
private:
  Solver<T> &m_solver;
  //  std::vector<task> m_tasks;
  Interval<T> schedule;
  std::vector<Interval<T>> the_tasks;
  std::vector<std::vector<Literal<T>>> disjunct;

  // helpers
  std::vector<unsigned> est_order;
  std::vector<unsigned> lct_order;
  std::vector<unsigned> theta_rank;
  ThetaTree TT;

  //  std::vector<std::vector<task>> explanation_tasks;
  std::vector<std::vector<Interval<T> *>> the_explanation_tasks;
  std::vector<T> explanation_lb;
  std::vector<T> explanation_ub;

  Reversible<size_t> num_explanations;

  std::vector<unsigned> pruned_tasks;
  std::vector<std::vector<unsigned>::reverse_iterator> omegas;
  std::vector<std::vector<unsigned>::iterator> bomegas;
  std::vector<T> relevant_starts;
  std::vector<T> bound_omegas;

  template <typename Iter>
  hint lowerBoundExplanation(Iter b, Iter e, const T lb);
  template <typename Iter>
  hint upperBoundExplanation(Iter b, Iter e, const T ub);

  //  bool falsified(const lit e);

  void propagateForward();
  void propagateBackward();

public:
  template <typename ItTask, typename ItVar>
  DisjunctiveEdgeFinding(Solver<T> &solver, Interval<T> &sched, ItTask beg_task,
                         ItTask end_task, ItVar beg_var);
  virtual ~DisjunctiveEdgeFinding();

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;

  void printLBExplanation(const hint ph);
  void printUBExplanation(const hint ph);
  void printTrivialExplanation(const Literal<T> l);

  //  static std::vector<int> task_map;
  //#ifdef DBG_EDGEFINDING
  //  int debug_flag{0};
  //#endif

  template <typename Iter>
  bool checklbpruning(const unsigned r, const T lb, const Iter b, const Iter e);

  template <typename Iter>
  bool checkubpruning(const unsigned r, const T lb, const Iter b, const Iter e);
};

template <typename T>
std::string DisjunctiveEdgeFinding<T>::prettyTask(const int i) const {
  std::stringstream ss;
  //  ss << "t" << m_tasks[i] << ": [" << est(i) << ".." << lct(i) << "] ("
  //     << minduration(i) << ")";
  ss << "t" << the_tasks[i].id() << ": [" << est(i) << ".." << lct(i) << "] ("
     << minduration(i) << ")";
  return ss.str();
}

// collect all the tasks in [b,e) whose lower bound is larger than or equal to
// lb
template <typename T>
template <typename Iter>
hint DisjunctiveEdgeFinding<T>::lowerBoundExplanation(Iter b, Iter e,
                                                      const T lb) {
  auto e_idx{num_explanations};
  if (the_explanation_tasks.size() <= e_idx) {
    the_explanation_tasks.resize(e_idx + 1);
    explanation_lb.resize(e_idx + 1);
    explanation_ub.resize(e_idx + 1);
  } else {
    the_explanation_tasks[e_idx].clear();
  }
  explanation_lb[e_idx] = lb;

  T ub{the_tasks[*b].getLatestEnd(m_solver)};
  for (auto x{b}; x != e; ++x) {
    if (the_tasks[*x].getEarliestStart(m_solver) >= lb) {

#ifdef DBG_EXPLEF
      std::cout << " expl: " the_tasks[*x].start.after(lb);
<<)
<< " and "
the_tasks[*x].end.before(the_tasks[*x].getLatestEnd(m_solver));
<< std::endl;
#endif

//      explanation_tasks[e_idx].push_back(i);
the_explanation_tasks[e_idx].push_back(&the_tasks[*x]);

ub = std::max(ub, the_tasks[*x].getLatestEnd(m_solver));
    }
#ifdef DBG_EXPLEF
    else {
      std::cout << " skip " << prettyTask(*x) << " b/c of the lower bound ("
                << lb << ")" << std::endl;
    }
#endif
  }

  assert(ub >= the_tasks[*b].getLatestEnd(m_solver));
  explanation_ub[e_idx] = ub;

  ++num_explanations;

  return static_cast<hint>(e_idx);
}

// collect all the tasks in [b,e) whose upper bound is lower than or equal to ub
template <typename T>
template <typename Iter>
hint DisjunctiveEdgeFinding<T>::upperBoundExplanation(Iter b, Iter e,
                                                      const T ub) {
  auto e_idx{num_explanations};
  if (the_explanation_tasks.size() <= e_idx) {
    //    explanation_tasks.resize(e_idx + 1);
    the_explanation_tasks.resize(e_idx + 1);
    explanation_lb.resize(e_idx + 1);
    explanation_ub.resize(e_idx + 1);
  } else {
    //    explanation_tasks[e_idx].clear();
    the_explanation_tasks[e_idx].clear();
  }
  explanation_ub[e_idx] = ub;

  T lb{the_tasks[*b].getEarliestStart(m_solver)};

  for (auto x{b}; x != e; ++x) {

    if (the_tasks[*x].getLatestEnd(m_solver) <= ub) {

#ifdef DBG_EXPLEF
      std::cout << " expl: " the_tasks[*x]->start.after(
          the_tasks[*x]->getEarliestStart(m_solver));
      << " and " << the_tasks[*x]->end.before(ub);
      << std::endl;
#endif
      the_explanation_tasks[e_idx].push_back(&the_tasks[*x]);
      lb = std::min(lb, the_tasks[*x].getEarliestStart(m_solver));
    }
  }

  assert(lb <= the_tasks[*b].getEarliestStart(m_solver));
  explanation_lb[e_idx] = lb;
  ++num_explanations; // = explanation_tasks.size();

  return static_cast<hint>(e_idx);
}

template <typename T> T DisjunctiveEdgeFinding<T>::est(const unsigned i) const {
  return the_tasks[i].getEarliestStart(m_solver);
}

template <typename T> T DisjunctiveEdgeFinding<T>::lst(const unsigned i) const {
  return the_tasks[i].getLatestStart(m_solver);
}

template <typename T> T DisjunctiveEdgeFinding<T>::ect(const unsigned i) const {
  return the_tasks[i].getEarliestEnd(m_solver);
}

template <typename T> T DisjunctiveEdgeFinding<T>::lct(const unsigned i) const {
  return the_tasks[i].getLatestEnd(m_solver);
}

template <typename T>
T DisjunctiveEdgeFinding<T>::minduration(const unsigned i) const {
  return the_tasks[i].minDuration();
}

template <typename T>
T DisjunctiveEdgeFinding<T>::maxduration(const unsigned i) const {
  return the_tasks[i].maxDuration();
}

template <typename T>
template <typename ItTask, typename ItVar>
DisjunctiveEdgeFinding<T>::DisjunctiveEdgeFinding(Solver<T> &solver,
                                                  Interval<T> &sched,
                                                  ItTask beg_task,
                                                  ItTask end_task,
                                                  ItVar beg_var)
    : m_solver(solver), TT(std::distance(beg_task, end_task)),
      num_explanations(0, &(m_solver.getEnv())) {

  schedule = sched,

  Constraint<T>::priority = Priority::Medium;

  for (auto j{beg_task}; j != end_task; ++j) {
    the_tasks.push_back(*j);
  }

  disjunct.resize(the_tasks.size());

  for (unsigned i = 0; i < the_tasks.size(); ++i) {
    lct_order.push_back(i);
    est_order.push_back(i);
    disjunct[i].resize(the_tasks.size());
  }

  auto ep{beg_var};
  for (auto ip{beg_task}; ip != end_task; ++ip) {
    for (auto jp{ip + 1}; jp != end_task; ++jp) {
      auto x{*ep};

      auto i{std::distance(beg_task, ip)};
      auto j{std::distance(beg_task, jp)};
      disjunct[i][j] = m_solver.boolean.getLiteral(true, x);
      disjunct[j][i] = m_solver.boolean.getLiteral(false, x);

      ++ep;
    }
  }

  theta_rank.resize(the_tasks.size(), 0);
}

template <typename T> DisjunctiveEdgeFinding<T>::~DisjunctiveEdgeFinding() {}

template <typename T> void DisjunctiveEdgeFinding<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (unsigned i{0}; i < the_tasks.size(); ++i) {
    m_solver.wake_me_on(lb<T>(the_tasks[i].getStart()), this->id());
    m_solver.wake_me_on(ub<T>(the_tasks[i].getEnd()), this->id());
  }
}

template <typename T>
template <typename Iter>
bool DisjunctiveEdgeFinding<T>::checklbpruning(const unsigned r, const T lb,
                                               const Iter b, const Iter e) {
  T ub{0};
  T duration{0};
  bool trivial{true};
  for (auto it{b}; it != e; ++it) {

    trivial &= est(r) + minduration(r) + minduration(*it) > lct(*it);

    if (lb <= est(*it)) {
      duration += minduration(*it);
      ub = std::max(ub, lct(*it));
    }
  }

  assert(est(r) >= lb);
  duration += minduration(r);

  return trivial or (duration > (ub - lb));
}

template <typename T>
template <typename Iter>
bool DisjunctiveEdgeFinding<T>::checkubpruning(const unsigned r, const T ub,
                                               const Iter b, const Iter e) {

  T lb{Constant::Infinity<T>};
  T duration{0};
  bool trivial{true};
  for (auto it{b}; it != e; ++it) {

    trivial &= est(*it) + minduration(r) + minduration(*it) > lct(r);

    if (ub >= lct(*it)) {
      duration += minduration(*it);
      lb = std::min(lb, est(*it));
    }
  }

  assert(lct(r) <= ub);
  duration += minduration(r);

  return trivial or (duration > (ub - lb));
}

template <typename T>
bool DisjunctiveEdgeFinding<T>::notify(const Literal<T>, const int) {
  return true;
}


template <typename T>
void DisjunctiveEdgeFinding<T>::printLBExplanation(const hint ph) {
  auto l{the_explanation_tasks[ph].back()};
  std::cout << "because t" << l->id() << " starts after " << explanation_lb[ph]
            << " and";
  T dur{l->minDuration()};
  for (auto ti : the_explanation_tasks[ph]) {
    if (ti != l) {
      std::cout << " t" << ti->id();
      std::cout.flush();
      dur += ti->minDuration();
      assert(ti->getEarliestStart(m_solver) >= explanation_lb[ph]);
      assert(ti->getLatestEnd(m_solver) <= explanation_ub[ph]);
    }
  }
  std::cout << " are all in [" << explanation_lb[ph] << ".."
            << explanation_ub[ph] << "] and together they have a duration of "
            << dur << " which is larger than " << explanation_ub[ph] << " - "
            << explanation_lb[ph] << " = "
            << (explanation_ub[ph] - explanation_lb[ph]) << std::endl;
  assert(l->getEarliestStart(m_solver) >= explanation_lb[ph]);
  assert(dur > (explanation_ub[ph] - explanation_lb[ph]));
}

template <typename T>
void DisjunctiveEdgeFinding<T>::printUBExplanation(const hint ph) {
  auto l{the_explanation_tasks[ph].back()};
  std::cout << " [" << explanation_lb[ph] << ".." << explanation_ub[ph]
            << "] because t" << l->id() << " ends before " << explanation_ub[ph]
            << " and";
  T dur{l->minDuration()};
  for (auto ti : the_explanation_tasks[ph]) {
    if (ti != l) {
      std::cout << " t" << ti->id();
      std::cout.flush();
      dur += ti->minDuration();

      assert(ti->getEarliestStart(m_solver) >= explanation_lb[ph]);
      assert(ti->getLatestEnd(m_solver) <= explanation_ub[ph]);
    }
  }
  std::cout << " are all in [" << explanation_lb[ph] << ".."
            << explanation_ub[ph] << "] and together they have a duration of "
            << dur << " which is larger than " << explanation_ub[ph] << " - "
            << explanation_lb[ph] << " = "
            << (explanation_ub[ph] - explanation_lb[ph]) << std::endl;
  assert(l->getLatestEnd(m_solver) <= explanation_ub[ph]);
  assert(dur > (explanation_ub[ph] - explanation_lb[ph]));
}

template <typename T>
void DisjunctiveEdgeFinding<T>::printTrivialExplanation(const Literal<T> l) {

  auto eij{m_solver.boolean.getEdge(l)};
  std::cout << " because x" << eij.from << " = "
            << m_solver.numeric.lower(eij.from) << " + "
            << "?"
            << " + "
            << "?"
            << " > " << m_solver.numeric.upper(eij.to) << " = x" << eij.to << std::endl;

  //  assert(m_schedule.lower(eij.from) + m_schedule.minDuration(TASK(eij.from))
  //  +
  //             m_schedule.minDuration(TASK(eij.to)) >
  //         m_schedule.upper(eij.to));
}

template <typename T> void DisjunctiveEdgeFinding<T>::propagate() {
  propagateForward();
  propagateBackward();
}

template <typename T> void DisjunctiveEdgeFinding<T>::propagateForward() {

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "propagate edge-finding forward(";
    std::cout << prettyTask(0);
    for (size_t i{1}; i < the_tasks.size(); ++i) {
        std::cout << " " << prettyTask(i);
    }
    std::cout << ") (i=" << m_solver.num_cons_propagations << ")\n";
  }
#endif

  hint ph{Constant::FactHint};

  std::sort(est_order.begin(), est_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return est(a) < est(b);
            });

  std::sort(lct_order.begin(), lct_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return lct(a) < lct(b);
            });

  TT.clear();

  for (unsigned i{0}; i < est_order.size(); ++i) {
    theta_rank[est_order[i]] = i;
  }

  // insert the tasks one by one in the theta tree by non-decreasing lct order
  for (auto ai{lct_order.begin()}; ai != lct_order.end(); ++ai) {
    auto a{*ai};
    TT.insert(theta_rank[a], est(a), minduration(a));

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "insert " << prettyTask(a) << " bound=" << TT.getBound()
                << std::endl;
    }
#endif

    // failure because of overload
    if (TT.getBound() > lct(a)) {
      auto h{lowerBoundExplanation(lct_order.begin(), ai + 1, TT.getEst())};
      throw Failure<T>({this, h});
    }
  }

  // global upper bound from the overload formula
  if (TT.getBound() > schedule.getEarliestEnd(m_solver)) {
    auto h{
        lowerBoundExplanation(lct_order.begin(), lct_order.end(), TT.getEst())};
    m_solver.set(schedule.end.after(-TT.getBound()), {this, h});
  }

  // Edge-finding
  pruned_tasks.clear();
  omegas.clear();
  relevant_starts.clear();
  bound_omegas.clear();

  // Edge-finding check on the tasks, one by on in the reverse insertion order
  for (auto ai{lct_order.rbegin()}; ai != (lct_order.rend() - 1); ++ai) {
    auto a{*ai};

    // deadline without the gray task
    auto deadline_omega{lct(*(ai + 1))};
    TT.paint_gray(theta_rank[a], a);

    // bound including the processing of the gray task
    auto ect_{TT.grayBound()};
    assert(TT.getBound() <= deadline_omega);
#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "gray " << prettyTask(a) << " lct=" << deadline_omega
                << ", bound=" << ect_ << std::endl;
    }
#endif

    // record the tasks that should be pruned (pruning online can mess up with
    // the ordering and is error-prone)
    while (ect_ > deadline_omega) {
      auto r{TT.getResponsible()};
      pruned_tasks.push_back(r);

      // the upper bound of the tasks' windows
      omegas.push_back(ai + 1);

      // the lower bound of the tasks' windows
      relevant_starts.push_back(TT.grayEst());

      // the new lower bound for r
      bound_omegas.push_back(ect_);

      // actually remove r
      TT.remove(theta_rank[r]);
      ect_ = TT.grayBound();
#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << prettyTask(r) << " must be after:\n";
        for (auto it{omegas.back()}; it != lct_order.rend(); ++it) {
          std::cout << " - " << prettyTask(*it) << std::endl;
        }
        std::cout << "rm " << prettyTask(r) << " bound=" << ect_ << std::endl;
      }
#endif
    }
  }

  // do the pruning
  while (not pruned_tasks.empty()) {
    auto r{pruned_tasks.back()};
    auto ai{omegas.back()};
    auto s{relevant_starts.back()};
    auto ect_omega{bound_omegas.back()};
    pruned_tasks.pop_back();
    omegas.pop_back();
    relevant_starts.pop_back();
    bound_omegas.pop_back();

    auto &tl{the_tasks[r]};
    ph = Constant::FactHint;
    for (auto j{ai}; j != lct_order.rend(); ++j) {

      if (not m_solver.boolean.satisfied(disjunct[r][*j])) {

#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          std::cout << "add precedence " << disjunct[r][*j] << "?" << std::endl;
        }
#endif

        // ej < si (ub(si) & lb(ej))
        if (m_solver.boolean.falsified(disjunct[*j][r])) {
          // the precedence is trivially implied because r >> *j
          // not clear if we should just let the edge constraints handle that

          // this is the edge that is satisfied and implies disjunct[r][*j]
          auto eij{~m_solver.boolean.getEdge(disjunct[*j][r])};
          auto idx_lb{m_solver.numeric.lastLitIndex(bound::lower, eij.from)};
          auto idx_ub{m_solver.numeric.lastLitIndex(bound::upper, eij.to)};

          // ei < sj (ub(ei) & lb(sj))
          m_solver.set(
              disjunct[r][*j],
              {this, -1 - static_cast<hint>(std::min(idx_lb, idx_ub))});

#ifdef DBG_EDGEFINDING
          if (DBG_EDGEFINDING) {
            printTrivialExplanation(disjunct[r][*j]);
          }
#endif
        } else {
          assert(checklbpruning(r, s, ai, lct_order.rend()));
          if (ph == Constant::FactHint) {
            ph = lowerBoundExplanation(ai, lct_order.rend(), s);
            //                        explanation_tasks[ph].push_back(l);
            the_explanation_tasks[ph].push_back(&tl);
#ifdef DBG_EDGEFINDING
            if (DBG_EDGEFINDING) {
              printLBExplanation(ph);
            }
#endif
          }
          m_solver.set(disjunct[r][*j], {this, ph});
        }
      }
    }

    if (ect(r) < ect_omega) {

      Literal<T> bc{tl.end.after(ect_omega)};

      //          BoundConstraint<T> bc_verif{LOWERBOUND(END(l)),
      //          -ect_omega};

      //          assert(bc == bc_verif);

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << "update bound " << bc << std::endl;
      }
#endif

      if (ph == Constant::FactHint) {
        ph = lowerBoundExplanation(ai, lct_order.rend(), s);
        //                explanation_tasks[ph].push_back(l);
        the_explanation_tasks[ph].push_back(&tl);

#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          printLBExplanation(ph);
        }
#endif
      }
      m_solver.set(bc, {this, ph});
    }
  }
}


template <typename T> void DisjunctiveEdgeFinding<T>::propagateBackward() {

  auto horizon(schedule.getLatestEnd(m_solver));

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "propagate edge-finding backward(";
    std::cout << prettyTask(0);
    for (size_t i{1}; i < the_tasks.size(); ++i) {
      std::cout << " " << prettyTask(i);
    }
    std::cout << ") (i=" << m_solver.num_cons_propagations << ")\n";
  }
#endif

  hint ph{Constant::FactHint};

  std::sort(est_order.begin(), est_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return est(a) < est(b);
            });

  std::sort(lct_order.begin(), lct_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return lct(a) < lct(b);
            });

  TT.clear();

  for (unsigned i{0}; i < lct_order.size(); ++i) {
    theta_rank[lct_order[lct_order.size() - i - 1]] = i;
  }

  for (auto ai{est_order.rbegin()}; ai != est_order.rend(); ++ai) {
    auto a{*ai};

    TT.insert(theta_rank[a], horizon - lct(a), minduration(a));

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "insert " << prettyTask(a)
                << " bound=" << (horizon - TT.getBound()) << std::endl;
    }
#endif

    if (TT.getBound() > (horizon - est(a))) {
      auto h{upperBoundExplanation(est_order.rbegin(), ai + 1,
                                   horizon - TT.getEst())};
      throw Failure<T>({this, h});
    }
  }

  pruned_tasks.clear();
  bomegas.clear();
  relevant_starts.clear();
  bound_omegas.clear();
  for (auto ai{est_order.begin()}; ai != est_order.end() - 1; ++ai) {
    auto a{*ai};
    auto deadline_omega{horizon - est(*(ai + 1))};
    TT.paint_gray(theta_rank[a], a);
    auto ect_{TT.grayBound()};
    assert(TT.getBound() <= deadline_omega);
#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "gray " << prettyTask(a)
                << " lct=" << (horizon - deadline_omega)
                << ", bound=" << (horizon - ect_) << std::endl;
    }
#endif
    while (ect_ > deadline_omega) {
      auto r{TT.getResponsible()};
      pruned_tasks.push_back(r);
      bomegas.push_back(ai);
      relevant_starts.push_back(horizon - TT.grayEst());
      bound_omegas.push_back(horizon - ect_);

      TT.remove(theta_rank[r]);
      ect_ = TT.grayBound();
#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << prettyTask(r) << " must be before:\n";
        for (auto it{bomegas.back() + 1}; it != est_order.end(); ++it) {
          std::cout << " - " << prettyTask(*it) << std::endl;
        }
        std::cout << "rm " << prettyTask(r) << " bound=" << ect_ << std::endl;
      }
#endif
    }
  }

  while (not pruned_tasks.empty()) {
    auto r{pruned_tasks.back()};
    auto ai{bomegas.back()};
    auto u{relevant_starts.back()};
    auto lst_omega{bound_omegas.back()};
    pruned_tasks.pop_back();
    bomegas.pop_back();
    relevant_starts.pop_back();
    bound_omegas.pop_back();

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "add precedences w.r.t. " << prettyTask(r) << std::endl;
    }
#endif

    //        auto f{m_tasks[r]};
    auto &fl{the_tasks[r]};
    ph = Constant::FactHint;
    for (auto j{est_order.rbegin()}; *j != *ai; ++j) {
      if (not m_solver.boolean.satisfied(disjunct[*j][r])) {
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          std::cout << "add precedence " << disjunct[*j][r] << "?" << std::endl;
        }
#endif
        // ej < si (ub(si) & lb(ej))
        if (m_solver.boolean.falsified(disjunct[r][*j])) {

          auto eij{~m_solver.boolean.getEdge(disjunct[r][*j])};
          auto idx_lb{m_solver.numeric.lastLitIndex(bound::lower, eij.from)};
          auto idx_ub{m_solver.numeric.lastLitIndex(bound::upper, eij.to)};

          // ei < sj (ub(ei) & lb(sj))
          m_solver.set(
              disjunct[*j][r],
              {this, -1 - static_cast<hint>(std::min(idx_lb, idx_ub))});

#ifdef DBG_EDGEFINDING
          if (DBG_EDGEFINDING) {
            printTrivialExplanation(disjunct[*j][r]);
          }
#endif
        } else {
          assert(checkubpruning(r, u, ai + 1, est_order.end()));
          if (ph == Constant::FactHint) {
            ph = upperBoundExplanation(ai + 1, est_order.end(), u);
            //                        explanation_tasks[ph].push_back(f);
            the_explanation_tasks[ph].push_back(&fl);
#ifdef DBG_EDGEFINDING
            if (DBG_EDGEFINDING) {
              printUBExplanation(ph);
            }
#endif
          }
          m_solver.set(disjunct[*j][r], {this, ph});
        }
      }
    }

    if (lst(r) < lst_omega) {

      Literal<T> bc{fl.start.before(lst_omega)};

      //          BoundConstraint<T> bc_verif{UPPERBOUND(START(f)),
      //          lst_omega};

      //          assert(bc == bc_verif);

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << "update bound " << bc << std::endl;
      }
#endif

      if (ph == Constant::FactHint) {
        ph = upperBoundExplanation(ai + 1, est_order.end(), u);
        //                explanation_tasks[ph].push_back(f);
        the_explanation_tasks[ph].push_back(&fl);

#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          printUBExplanation(ph);
        }
#endif
      }
      m_solver.set(bc, {this, ph});
    }
  }
}

template <typename T> int DisjunctiveEdgeFinding<T>::getType() const {
  return EDGEFINDINGEXPL;
}

template <typename T>
void DisjunctiveEdgeFinding<T>::xplain(const Literal<T> l, const hint h,
                                       std::vector<Literal<T>> &Cl) {

  //  if (static_cast<size_t>(h) >= explanation_tasks.size()) {
  //    std::cout << h << " / " << explanation_tasks.size() << std::endl;
  //    exit(1);
  //  }

  if (l == Solver<T>::Contradiction) {
#ifdef DBG_EXPLEF
    std::cout << "explain failure from edge-finding: overload on interval ["
              << explanation_lb[h] << ".." << explanation_ub[h] << "]\n";
    T duration{0};
#endif

    // failure case, everything is in "tasks"
    for (auto ti : the_explanation_tasks[h]) {

      Cl.push_back(ti->start.after(explanation_lb[h]));
      Cl.push_back(ti->end.before(explanation_ub[h]));

#ifdef DBG_EXPLEF
      duration += ti->minDuration();
      std::cout << ti->start.after(explanation_lb[h]) << " & "
                << ti->end.before(explanation_ub[h]) << " ("
                << ti->minDuration() << ")\n";
#endif
    }

#ifdef DBG_EXPLEF
    assert(duration > (explanation_ub[h] - explanation_lb[h]));
#endif
  } else if (h < 0) {

    //#ifdef DBG_EXPLEF
    std::cout << "explain " << m_solver.pretty(l)
              << " (binary disjunctive reasoning) TODO!!\n";
    //#endif

    //    auto exy{m_solver.boolean.getEdge(l)};
    //    auto bidx{-1 - h};

    exit(1);

    //    auto lbl = -1 - h;
    //    auto lbc = m_schedule.getBound(lbl);
    //
    //    //        std::cout << "edge: " << exy << std::endl;
    //    //        std::cout << "lb (" << lbl << "): " << lbc << " / " <<
    //    //        m_schedule.prettyLiteral(BOUND(lbl)) << std::endl;
    //    //
    //
    //    //        std::cout << lbc.distance << " + " <<
    //    //        m_schedule.minDuration(TASK(exy.to)) << " + " <<
    //    //        m_schedule.minDuration(TASK(exy.from)) << std::endl;
    //
    //    //    BoundConstraint<T> ubc{UPPERBOUND(exy.to),
    //    //                           m_schedule.minDuration(TASK(exy.to)) +
    //    // m_schedule.minDuration(TASK(exy.from)) -
    //    //                               lbc.distance - Gap<T>::epsilon()};
    //
    //    BoundConstraint<T> ubc{UPPERBOUND(exy.to),
    //                           m_schedule.getTask(exy.to).minDuration() +
    //                               // m_schedule.minDuration(TASK(exy.to))
    //                               //                             +
    //                               m_schedule.getTask(exy.from).minDuration()
    //                               -
    //                               // m_schedule.minDuration(TASK(exy.from))
    //                               //                                 -
    //                               lbc.distance - Gap<T>::epsilon()};
    //
    //    auto ubl{m_schedule.getImplicant(ubc)};
    //
    //#ifdef DBG_EXPLEF
    //    std::cout << m_schedule.prettyLiteral(BOUND(lbl)) << " and "
    //              << m_schedule.prettyLiteral(BOUND(ubl)) << std::endl;
    //#endif
    //
    //    Cl.push_back(BOUND(lbl));
    //    Cl.push_back(BOUND(ubl));
    //
    //    //        exit(1);
    //    //#endif

  } else {

#ifdef DBG_EXPLEF
    std::cout << "explain " << m_solver.pretty(l)
              << " (overload reasoning on interval [" << explanation_lb[h]
              << ".." << explanation_ub[h] << "])\n";
    T duration{0};
#endif

    // edge case, "tasks" give the Omega
    auto n{the_explanation_tasks[h].size() - 1};
    for (size_t i{0}; i < n; ++i) {
      auto ti{the_explanation_tasks[h][i]};

      Cl.push_back(ti->start.after(explanation_lb[h]));
      Cl.push_back(ti->end.before(explanation_ub[h]));

#ifdef DBG_EXPLEF
      //          duration += m_schedule.minDuration(t);
      duration += ti->minDuration();
      std::cout << ti->start.after(explanation_lb[h]) << " & "
                << ti->end.before(explanation_ub[h]) << " ("
                << ti->minDuration();
      //          m_schedule.minDuration(t)
      << ")\n";
#endif
    }

    // the last task is the one at the tip of the edges
    auto t{the_explanation_tasks[h].back()};

    if (l.isNumeric()) {
      // explain the bound change
      if (l.sign() == bound::lower) {
        Cl.push_back(t->start.after(explanation_lb[h]));
      } else {
        Cl.push_back(t->end.before(explanation_ub[h]));
      }
    } else {
      // explain the edges
      auto lc{m_solver.boolean.getEdge(l)};
      if (t->getStart() == lc.from or t->getEnd() == lc.from) {
        Cl.push_back(t->start.after(explanation_lb[h]));
      } else {
        Cl.push_back(t->end.before(explanation_ub[h]));
      }
    }

#ifdef DBG_EXPLEF
    duration += t->minDuration(); // m_schedule.minDuration(i);
    if (t.getStart() == lc.from or t.getEnd() == lc.from) {
      std::cout << t->start.after(explanation_lb[h]);
    } else {
      std::cout << t->end.before(explanation_ub[h]);
    }
    std::cout << " (" << t->minDuration() << ")\n";

    std::cout << duration << " > " << explanation_ub[h] << " - "
              << explanation_lb[h] << std::endl;
    assert(duration > (explanation_ub[h] - explanation_lb[h]));
#endif
  }
}

template <typename T>
std::ostream &DisjunctiveEdgeFinding<T>::display(std::ostream &os) const {
  os << "Disjunctive Edge-Finding";

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
std::ostream &DisjunctiveEdgeFinding<T>::print_reason(std::ostream &os,
                                                      const hint) const {
  os << "edge-finding";
  return os;
}

} // namespace tempo

#endif
