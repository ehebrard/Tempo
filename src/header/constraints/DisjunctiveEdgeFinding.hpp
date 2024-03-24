#ifndef TEMPO_DISJUNCTIVEEDGEFINDING_HPP
#define TEMPO_DISJUNCTIVEEDGEFINDING_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Scheduler.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"
#include "util/ThetaTree.hpp"


namespace tempo {

template <typename T> class DisjunctiveEdgeFinding : public Constraint {
private:
  Scheduler<T> &m_schedule;
  std::vector<task> m_tasks;
  std::vector<std::vector<lit>> disjunct;

  // helpers
  std::vector<unsigned> est_order;
  std::vector<unsigned> lct_order;
  std::vector<unsigned> theta_rank;
  ThetaTree TT;

  std::vector<std::vector<task>> nexplanation_tasks;
  std::vector<T> nexplanation_lb;
  std::vector<T> nexplanation_ub;

  Reversible<size_t> num_explanations;

  template <typename Iter>
  hint newLowerBoundExplanation(Iter b, Iter e, const T lb);
  template <typename Iter>
  hint newUpperBoundExplanation(Iter b, Iter e, const T ub);

public:
  template <typename ItTask, typename ItVar>
  DisjunctiveEdgeFinding(Scheduler<T> &scheduler, const ItTask beg_task,
                         const ItTask end_task, const ItVar beg_var,
                         const ItVar end_var);
  virtual ~DisjunctiveEdgeFinding();

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;

  bool notify_bound(const int lit, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  static std::vector<int> task_map;
  //#ifdef DBG_EDGEFINDING
  //  int debug_flag{0};
  //#endif
};



// collect all the tasks in [b,e) whose lower bound is larger than or equal to
// lb
template <typename T>
template <typename Iter>
hint DisjunctiveEdgeFinding<T>::newLowerBoundExplanation(Iter b, Iter e,
                                                         const T lb) {
  auto e_idx{num_explanations};
  if (nexplanation_tasks.size() <= e_idx) {
    nexplanation_tasks.resize(e_idx + 1);
    nexplanation_lb.resize(e_idx + 1);
    nexplanation_ub.resize(e_idx + 1);
  } else {
    nexplanation_tasks[e_idx].clear();
  }
  nexplanation_lb[e_idx] = lb;

  T ub{m_schedule.upper(END(m_tasks[*b]))};
  for (auto x{b}; x != e; ++x) {

    auto i{m_tasks[*x]};

    if (m_schedule.lower(START(i)) >= lb) {

#ifdef DBG_EXPLEF
      std::cout << " expl: "
                << m_schedule.prettyLiteral(
                       BOUND(m_schedule.getBoundIndex(LOWERBOUND(START(i)))))
                << " and "
                << m_schedule.prettyLiteral(
                       BOUND(m_schedule.getBoundIndex(UPPERBOUND(END(i)))))
                << std::endl;
#endif

      nexplanation_tasks[e_idx].push_back(i);
      ub = std::max(ub, m_schedule.upper(END(i)));
    }
  }
  assert(ub >= m_schedule.upper(END(m_tasks[*b])));
  nexplanation_ub[e_idx] = ub;

  ++num_explanations;

  return static_cast<hint>(e_idx);
}

// collect all the tasks in [b,e) whose upper bound is lower than or equal to ub
template <typename T>
template <typename Iter>
hint DisjunctiveEdgeFinding<T>::newUpperBoundExplanation(Iter b, Iter e,
                                                         const T ub) {
  auto e_idx{num_explanations};
  if (nexplanation_tasks.size() <= e_idx) {
    nexplanation_tasks.resize(e_idx + 1);
    nexplanation_lb.resize(e_idx + 1);
    nexplanation_ub.resize(e_idx + 1);
  } else {
    nexplanation_tasks[e_idx].clear();
  }
  nexplanation_ub[e_idx] = ub;

  T lb{m_schedule.lower(START(m_tasks[*b]))};
  for (auto x{b}; x != e; ++x) {

    auto i{m_tasks[*x]};

    if (m_schedule.lower(END(i)) <= ub) {

#ifdef DBG_EXPLEF
      std::cout << " expl: "
                << m_schedule.prettyLiteral(
                       BOUND(m_schedule.getBoundIndex(LOWERBOUND(START(i)))))
                << " and "
                << m_schedule.prettyLiteral(
                       BOUND(m_schedule.getBoundIndex(UPPERBOUND(END(i)))))
                << std::endl;
#endif
      nexplanation_tasks[e_idx].push_back(i);
      lb = std::min(lb, m_schedule.lower(START(i)));
    }
  }

  assert(lb <= m_schedule.lower(START(m_tasks[*b])));
  nexplanation_lb[e_idx] = lb;

  ++num_explanations; // = nexplanation_tasks.size();

  return static_cast<hint>(e_idx);
}


template <typename T> T DisjunctiveEdgeFinding<T>::est(const unsigned i) const {
  return m_schedule.lower(START(m_tasks[i]));
}

template <typename T> T DisjunctiveEdgeFinding<T>::lst(const unsigned i) const {
  return m_schedule.upper(START(m_tasks[i]));
}

template <typename T> T DisjunctiveEdgeFinding<T>::ect(const unsigned i) const {
  return m_schedule.lower(END(m_tasks[i]));
}

template <typename T> T DisjunctiveEdgeFinding<T>::lct(const unsigned i) const {
  return m_schedule.upper(END(m_tasks[i]));
}

template <typename T>
T DisjunctiveEdgeFinding<T>::minduration(const unsigned i) const {
  return m_schedule.minDuration(m_tasks[i]);
}

template <typename T>
T DisjunctiveEdgeFinding<T>::maxduration(const unsigned i) const {
  return m_schedule.maxDuration(m_tasks[i]);
}

template <typename T>
template <typename ItTask, typename ItVar>
DisjunctiveEdgeFinding<T>::DisjunctiveEdgeFinding(Scheduler<T> &scheduler,
                                                  const ItTask beg_task,
                                                  const ItTask end_task,
                                                  const ItVar beg_var,
                                                  const ItVar end_var)
    : m_schedule(scheduler), TT(std::distance(beg_task, end_task)),
      num_explanations(0, &(m_schedule.getEnv())) {

  task_map.resize(m_schedule.numTask());

  // get all tasks with non-zero duration
  auto i{0};
  for (auto j{beg_task}; j != end_task; ++j) {

    task t{*j};
    task_map[t] = i++;
    m_tasks.push_back(t);

  }

  disjunct.resize(m_tasks.size());

  for (unsigned i = 0; i < m_tasks.size(); ++i) {
    lct_order.push_back(i);
    est_order.push_back(i);
    disjunct[i].resize(m_tasks.size());
  }

  for (auto v{beg_var}; v != end_var; ++v) {
    auto ep{m_schedule.getEdge(POS(*v))};
    auto pf{TASK(ep.from)};
    auto pt{TASK(ep.to)};

    auto en{m_schedule.getEdge(NEG(*v))};
    auto nf{TASK(en.from)};
    auto nt{TASK(en.to)};

    disjunct[task_map[pf]][task_map[pt]] = POS(*v);
    disjunct[task_map[nf]][task_map[nt]] = NEG(*v);
  }

  theta_rank.resize(m_tasks.size(), 0);

}

template <typename T> DisjunctiveEdgeFinding<T>::~DisjunctiveEdgeFinding() {}

template <typename T> void DisjunctiveEdgeFinding<T>::post(const int idx) {

  cons_id = idx;
  idempotent = true;

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (unsigned i{0}; i < m_tasks.size(); ++i) {
    m_schedule.wake_me_on_event(LOWERBOUND(START(m_tasks[i])), cons_id);
    m_schedule.wake_me_on_event(UPPERBOUND(END(m_tasks[i])), cons_id);
  }
}

template <typename T>
bool DisjunctiveEdgeFinding<T>::notify_bound(const lit, const int) {
  return true;
}

template <typename T> void DisjunctiveEdgeFinding<T>::propagate() {

//  bool pruning{false};
    hint ph{NoHint};

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "propagate edge-finding(";
    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << " t" << m_tasks[i] << " [" << est(i) << ".." << lct(i)
                << "]";
    }
    std::cout << ") (i=" << m_schedule.num_cons_propagations << ")\n";
  }
#ifdef DBG_THETA
  size_t debug_limit{0};
#endif
#endif

  std::sort(est_order.begin(), est_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return est(a) < est(b);
            });

  std::sort(lct_order.begin(), lct_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return lct(a) < lct(b);
            });

  auto horizon(m_schedule.upper(HORIZON));

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {

    for (auto i : est_order)
      std::cout << i << " ";
    std::cout << std::endl;

    for (auto i : lct_order)
      std::cout << i << " ";
    std::cout << std::endl;

    for (auto i : lct_order)
      std::cout << m_tasks[i] << " ";
    std::cout << std::endl;

    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "] (" << minduration(i) << ")\n";
    }
    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "] (" << minduration(i) << ")\n";
    }
    std::cout << std::endl << "forward:\n";
  }
#endif

  TT.clear();

  for (unsigned i{0}; i < est_order.size(); ++i) {
    theta_rank[est_order[i]] = i;
  }

  for (auto ai{lct_order.begin()}; ai != lct_order.end(); ++ai) {
    auto a{*ai};

    TT.insert(theta_rank[a], est(a), minduration(a));

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "add t" << m_tasks[a] << ": [" << est(a) << ".."
                << minduration(a) << ".." << lct(a) << "] : " << TT.getBound()
                << std::endl;
#ifdef DBG_THETA
      std::cout << TT << std::endl;
#endif
    }
#endif

    if (TT.getBound() > lct(a)) {
#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING)
        std::cout << "[efpruning] failure b/c overload check\n";
#endif
      auto h{newLowerBoundExplanation(lct_order.begin(), ai + 1, TT.getEst())};

#ifdef NO_REV
      throw Failure({this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
        throw Failure({this, h});
#endif
        
    }
  }

  if (TT.getBound() > m_schedule.lower(HORIZON)) {

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING)
      std::cout << "[efpruning] Cmax lower bound (" << -TT.getBound()
                << ") b/c overload check\n";
#endif

    auto h{newLowerBoundExplanation(lct_order.begin(), lct_order.end(),
                                    TT.getEst())};

#ifdef NO_REV
    m_schedule.set({LOWERBOUND(HORIZON), -TT.getBound()},
                   {this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
      m_schedule.set({LOWERBOUND(HORIZON), -TT.getBound()},
                     {this, h});
#endif
      
  }

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "] (" << minduration(i) << ")\n";
    }
    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "] (" << minduration(i) << ")\n";
    }
    std::cout << " forward edges:\n";
  }
#endif
    
  for (auto ai{lct_order.rbegin()}; ai != (lct_order.rend() - 1); ++ai) {
    auto a{*ai};
    auto deadline_omega{lct(*(ai + 1))};
    TT.paint_gray(theta_rank[a], a);

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "rm t" << m_tasks[a] << ": [" << est(a) << ".."
                << minduration(a) << ".." << lct(a) << "] : " << TT.getBound()
                << " (" << ect(a) << ")" << std::endl;
#ifdef DBG_THETA
      std::cout << TT << std::endl;
#endif
    }
#endif

    auto ect_{TT.grayBound()};

    if (TT.getBound() > deadline_omega) {

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING)
        std::cout << "[efpruning] failure b/c overload check\n";
#endif

      auto h{newLowerBoundExplanation(ai, lct_order.rend(), TT.getEst())};

#ifdef NO_REV
        throw Failure({this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
        throw Failure({this, h});
#endif
        
    }

#ifdef DBG_THETA
    debug_limit = m_tasks.size();
#endif

    while (ect_ > deadline_omega) {

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << "[efpruning] tasks";
        for (auto j{lct_order.begin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t" << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << " can't end before " << ect_ << ", but only t"
                  << m_tasks[r] << " can go past that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r) << ".." << lct(r)
                  << "] >= " << ect_ << " (" << TT.getEst() << "/"
                  << TT.grayEst() << ")" << std::endl;
      }
#endif

      auto l{m_tasks[r]};

//      hint ph{NoHint};
//      pruning = false;
        ph = NoHint;
      for (auto j{lct_order.begin()}; *j != a; ++j) {

        if (not m_schedule.satisfied(disjunct[r][*j])) {

          if (ph == NoHint) {
            ph =
                newLowerBoundExplanation(ai + 1, lct_order.rend(), TT.grayEst());
            nexplanation_tasks[ph].push_back(l);

//            pruning = true;
          }

#ifdef DBG_EDGEFINDING
          if (DBG_EDGEFINDING) {
            std::cout << "edge prec "
                      << m_schedule.prettyLiteral(EDGE(disjunct[r][*j]))
                      << std::endl;
          }
#endif

#ifdef NO_REV
            m_schedule.set(
                disjunct[r][*j],
                {this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
          m_schedule.set(
              disjunct[r][*j],
                         {this, ph});
#endif
        }
      }

      if (ect(r) < ect_) {

        BoundConstraint<T> bc{LOWERBOUND(END(l)), -ect_};
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          std::cout << "lower bound " << bc << std::endl;
        }
#endif
          
          if (ph == NoHint) {//} pruning) {
            ph =
                newLowerBoundExplanation(ai + 1, lct_order.rend(), TT.grayEst());
            nexplanation_tasks[ph].push_back(l);
//            pruning = true;
          }

#ifdef NO_REV
          m_schedule.set(
              bc, {this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
          m_schedule.set(
                         bc, {this, ph});
#endif
      }

//      assert((TT.grayEst() < est(r)) or not pruning);
        assert((TT.grayEst() < est(r)) or ph == NoHint);

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING and ph != NoHint) {
        std::cout << " because";
        for (auto l : explanation_tasks.back()) {
          std::cout << " t" << l;
        }
        std::cout << std::endl;
      }
#endif

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DBG_THETA
      std::cout << TT << std::endl;
      if (--debug_limit <= 0)
        exit(1);
      std::cout << debug_limit << std::endl;
#endif
    }
  }

  std::sort(est_order.begin(), est_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return est(a) < est(b);
            });

  std::sort(lct_order.begin(), lct_order.end(),
            [&](const unsigned a, const unsigned b) -> bool {
              return lct(a) < lct(b);
            });

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "] (" << minduration(i) << ")\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "] (" << minduration(i) << ")\n";
    }

    std::cout << "backward:\n";
  }
#endif

  TT.clear();
  for (unsigned i{0}; i < lct_order.size(); ++i) {
    theta_rank[lct_order[lct_order.size() - i - 1]] = i;
  }

  for (auto ai{est_order.rbegin()}; ai != est_order.rend(); ++ai) {
    auto a{*ai};

    TT.insert(theta_rank[a], horizon - lct(a), minduration(a));

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "add t" << m_tasks[a] << ": [" << (horizon - lct(a)) << ".."
                << minduration(a) << ".." << (horizon - est(a))
                << "] : " << TT.getBound() << " (" << (horizon - lst(a)) << ")"
                << std::endl;
#ifdef DBG_THETA
      std::cout << TT << std::endl;
#endif
    }
#endif
      
      if (TT.getBound() > (horizon - est(a))) {
  #ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING)
          std::cout << "[efpruning] failure b/c overload check\n";
  #endif
        auto h{newUpperBoundExplanation(est_order.rbegin(), ai + 1,
                                        horizon - TT.getEst())};

#ifdef NO_REV
          throw Failure({this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
          throw Failure({this, h});
#endif
      }
      
  }

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << m_schedule << "\nST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "] (" << minduration(i) << ")\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "] (" << minduration(i) << ")\n";
    }

    std::cout << " backward edges:\n";
  }
#endif

  for (auto ai{est_order.begin()}; ai != est_order.end() - 1; ++ai) {
    auto a{*ai};

    auto deadline_omega{horizon - est(*(ai + 1))};

    TT.paint_gray(theta_rank[a], a);

#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
      std::cout << "rm t" << m_tasks[a] << ": [" << (horizon - lct(a)) << ".."
                << minduration(a) << ".." << (horizon - est(a))
                << "] : " << TT.getBound() << " (" << (horizon - lst(a)) << ")"
                << std::endl;
#ifdef DBG_THETA
      std::cout << TT << std::endl;
#endif
    }
#endif

    auto ect_{TT.grayBound()};

    if (deadline_omega < TT.getBound()) {

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING)
        std::cout << "[efpruning] failure b/c overload check\n";
#endif

      auto h{
          newUpperBoundExplanation(ai, est_order.end(), horizon - TT.getEst())};

#ifdef NO_REV
        throw Failure({this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
        throw Failure({this, h});
#endif
    }

#ifdef DBG_THETA
    debug_limit = m_tasks.size();
#endif

    while (ect_ > deadline_omega) {

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << "[efpruning] tasks";
        for (auto j{est_order.rbegin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t" << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING) {
        std::cout << " can't start after " << (horizon - ect_) << ", but only t"
                  << m_tasks[r] << " can start at that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r) << ".." << lct(r)
                  << "] <= " << (horizon - ect_) << std::endl;
      }
#endif

      auto f{m_tasks[r]};

//      hint ph{NoHint};
//      pruning = false;
        ph = NoHint;
      for (auto j{est_order.rbegin()}; *j != a; ++j) {

        if (not m_schedule.satisfied(disjunct[*j][r])) {
          ++m_schedule.num_cons_prunings;
            if (ph == NoHint) { //not pruning) {
            ph = newUpperBoundExplanation(ai + 1, est_order.end(),
                                          horizon - TT.grayEst());
            nexplanation_tasks[ph].push_back(f);
//            pruning = true;
          }
#ifdef DBG_EDGEFINDING
          if (DBG_EDGEFINDING) {
            std::cout << "edge prec "
                      << m_schedule.prettyLiteral(EDGE(disjunct[r][*j]))
                      << std::endl;
          }
#endif

#ifdef NO_REV
            m_schedule.set(
                disjunct[*j][r],
                {this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
          m_schedule.set(
              disjunct[*j][r],
                         {this, ph});
#endif
        }
      }

      if (lst(r) > (horizon - ect_)) {

          if (ph == NoHint) { //not pruning) {
            ph = newUpperBoundExplanation(ai + 1, est_order.end(),
                                          horizon - TT.grayEst());
            nexplanation_tasks[ph].push_back(f);
//            pruning = true;
          }
          
        BoundConstraint<T> bc{UPPERBOUND(START(f)), (horizon - ect_)};
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
          std::cout << "upper bound " << bc << std::endl;
        }
#endif

#ifdef NO_REV
          m_schedule.set({UPPERBOUND(START(f)), (horizon - ect_)},
            {this, static_cast<hint>(explanation_tasks.size() - 1)});
#else
        m_schedule.set({UPPERBOUND(START(f)), (horizon - ect_)},
                       {this, ph});
#endif
      }

#ifdef DBG_EDGEFINDING
      if (DBG_EDGEFINDING and ph != NoHint) {
        std::cout << " because";
        for (auto l : explanation_tasks.back()) {
          std::cout << " t" << l;
        }
        std::cout << std::endl;
      }
#endif

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DBG_THETA
      std::cout << TT << std::endl;
      if (--debug_limit <= 0)
        exit(1);
      std::cout << debug_limit << std::endl;
#endif
    }
  }

#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "] (" << minduration(i) << ")\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "] (" << minduration(i) << ")\n";
    }
  }
#endif
}

template <typename T> int DisjunctiveEdgeFinding<T>::getType() const {
  return EDGEFINDINGEXPL;
}

template <typename T>
void DisjunctiveEdgeFinding<T>::xplain(const lit l, const hint h,
                                       std::vector<lit> &Cl) {
 
    if(static_cast<size_t>(h) >= nexplanation_tasks.size()) {
        std::cout << h << " / " << nexplanation_tasks.size() << std::endl;
        exit(1);
    }
    
    if (l == NoLit) {
#ifdef DBG_EXPLEF
        std::cout << "explain failure from edge-finding: overload on interval ["
        << nexplanation_lb[h] << ".." << nexplanation_ub[h] << "]\n";
        T duration{0};
#endif
        
        for (auto i : nexplanation_tasks[h]) {
            BoundConstraint<T> lb{LOWERBOUND(START(i)), -nexplanation_lb[h]};
            auto ll{m_schedule.getImplicant(lb)};
            Cl.push_back(BOUND(ll));
            
            BoundConstraint<T> ub{UPPERBOUND(END(i)), nexplanation_ub[h]};
            auto ul{m_schedule.getImplicant(ub)};
            Cl.push_back(BOUND(ul));
            
#ifdef DBG_EXPLEF
            duration += m_schedule.minDuration(i);
            std::cout << m_schedule.prettyLiteral(BOUND(ll)) << " & "
            << m_schedule.prettyLiteral(BOUND(ul)) << " ("
            << m_schedule.minDuration(i) << ")\n";
#endif
        }
        
#ifdef DBG_EXPLEF
        assert(duration > (nexplanation_ub[h] - nexplanation_lb[h]));
#endif
    } else {
        
#ifdef DBG_EXPLEF
        std::cout << "explain " << m_schedule.prettyLiteral(l) << " (overload reasoning on interval ["
        << nexplanation_lb[h] << ".." << nexplanation_ub[h] << "])\n";
        T duration{0};
#endif
        
        // failure case, everything is in "tasks"
        auto n{nexplanation_tasks[h].size() - 1};
        for (size_t i{0}; i < n; ++i) {
            task t{nexplanation_tasks[h][i]};
            
            BoundConstraint<T> lb{LOWERBOUND(START(t)), -nexplanation_lb[h]};
            auto ll{m_schedule.getImplicant(lb)};
            Cl.push_back(BOUND(ll));
            
            BoundConstraint<T> ub{UPPERBOUND(END(t)), nexplanation_ub[h]};
            auto ul{m_schedule.getImplicant(ub)};
            Cl.push_back(BOUND(ul));
            
#ifdef DBG_EXPLEF
            duration += m_schedule.minDuration(t);
            std::cout << m_schedule.prettyLiteral(BOUND(ll)) << " & "
            << m_schedule.prettyLiteral(BOUND(ul)) << " ("
            << m_schedule.minDuration(t) << ")\n";
#endif
        }
        
        //      bool bug{false};
        auto i{nexplanation_tasks[h].back()};
        lit p;
        if (LTYPE(l) == BOUND_LIT) {
            auto lc{m_schedule.getBound(FROM_GEN(l))};
            BoundConstraint<T> pc;
            if (SIGN(lc.l) == LOWER) {
                pc = {LOWERBOUND(START(i)), -nexplanation_lb[h]};
            } else {
                pc = {UPPERBOUND(END(i)), nexplanation_ub[h]};
            }
            p = m_schedule.getImplicant(pc);
            Cl.push_back(BOUND(p));
        } else {
            auto lc{m_schedule.getEdge(FROM_GEN(l))};
            BoundConstraint<T> pc;
            if (TASK(lc.from) == i) {
                pc = {LOWERBOUND(START(i)), -nexplanation_lb[h]};
            } else {
                pc = {UPPERBOUND(END(i)), nexplanation_ub[h]};
            }
            p = m_schedule.getImplicant(pc);
            Cl.push_back(BOUND(p));
        }
        
#ifdef DBG_EXPLEF
        duration += m_schedule.minDuration(i);
        std::cout << m_schedule.prettyLiteral(BOUND(p)) << " ("
        << m_schedule.minDuration(i) << ")\n";
        
        std::cout << duration << " > " << nexplanation_ub[h] << " - " << nexplanation_lb[h] << std::endl;
        assert(duration > (nexplanation_ub[h] - nexplanation_lb[h]));
#endif
        
    }
     
}

template <typename T>
std::ostream &DisjunctiveEdgeFinding<T>::display(std::ostream &os) const {
  os << "Disjunctive Edge-Finding";

#ifdef DBG_EDGEFINDING
  os << "[" << cons_id << "]";
#endif

  os << "(";
  for (auto t : m_tasks) {
    std::cout << " t" << t;
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &DisjunctiveEdgeFinding<T>::print_reason(std::ostream &os,
                                                      const hint) const {
  //  display(os);
  os << "edge-finding";
  //
  //  if (not explanations[h].empty()) {
  //
  //    auto l{explanations[h].begin()};
  //    m_schedule.displayLiteral(os, *l);
  //    ++l;
  //    while (l != explanations[h].end()) {
  //      os << ", ";
  //      m_schedule.displayLiteral(os, *l);
  //      ++l;
  //    }
  //  }
  //
  //  os << ")";
  return os;
}

template <typename T> std::vector<int> DisjunctiveEdgeFinding<T>::task_map;

} // namespace tempo

#endif
