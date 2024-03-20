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

  std::vector<lit> cur_explanation;
  std::vector<std::vector<lit>> explanations;

  // helpers
  std::vector<unsigned> est_order;
  std::vector<unsigned> lct_order;
  std::vector<unsigned> theta_rank;
  ThetaTree TT;

  template <typename Iter> void boundExplanation(Iter b, Iter e);

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

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  static std::vector<int> task_map;
#ifdef DEBUG_CONSTRAINT
  int debug_flag{0};
#endif
};

template <typename T>
template <typename Iter>
void DisjunctiveEdgeFinding<T>::boundExplanation(Iter b, Iter e) {

  assert(cur_explanation.empty());

  for (auto x{b}; x != e; ++x) {

    auto i{m_tasks[*x]};

    //      std::cout << "expl t" << i << ": "
    //      <<
    //      m_schedule.prettyLiteral(BOUND(m_schedule.getBoundIndex(LOWERBOUND(START(i)))))
    //      << " and "
    //      <<
    //      m_schedule.prettyLiteral(BOUND(m_schedule.getBoundIndex(UPPERBOUND(START(i)))))
    //      << " (" <<
    //      m_schedule.stamp(m_schedule.getBoundIndex(LOWERBOUND(START(i)))) <<
    //      ")" << std::endl;

    cur_explanation.push_back(
        BOUND(m_schedule.getBoundIndex(LOWERBOUND(START(i)))));
    cur_explanation.push_back(
        BOUND(m_schedule.getBoundIndex(UPPERBOUND(END(i)))));
  }

  explanations.push_back(cur_explanation);
  //  ++m_schedule.num_cons_prunings;
  cur_explanation.clear();
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
    : m_schedule(scheduler), TT(std::distance(beg_task, end_task)) {

  task_map.resize(m_schedule.numTask());

  // get all tasks with non-zero duration
  auto i{0};
  for (auto j{beg_task}; j != end_task; ++j) {

    task t{*j};
    task_map[t] = i++;
    //
    //      std::cout << "11\n";
    //
    //      std::cout <<

    if (m_schedule.maxDuration(t) > 0) {
      m_tasks.push_back(t);
    }
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

    //      std::cout << ep << std::endl << en << std::endl;
    //
    //    std::cout << *v << ": tasks t" << pf << " and t" << pt << " (" << nf
    //    << "/" << nt << "); indices " << task_map[pf]
    //              << " and " << task_map[pt] << " (" << task_map[nf] << "/" <<
    //              task_map[nt] <<"); lits " << POS(*v) << " and " << NEG(*v)
    //              << "; cons "
    //              << m_schedule.prettyLiteral(EDGE(POS(*v))) << " and " <<
    //              m_schedule.prettyLiteral(EDGE(NEG(*v))) << std::endl;

    disjunct[task_map[pf]][task_map[pt]] = POS(*v);
    disjunct[task_map[nf]][task_map[nt]] = NEG(*v);
  }

  theta_rank.resize(m_tasks.size(), 0);

  //        exit(1);

#ifdef DEBUG_CONSTRAINT
  debug_flag = 2;
#endif
}

template <typename T> DisjunctiveEdgeFinding<T>::~DisjunctiveEdgeFinding() {}

template <typename T> void DisjunctiveEdgeFinding<T>::post(const int idx) {

  cons_id = idx;
  idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
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

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << "propagate edge-finding(";
    for (size_t i{0}; i < m_tasks.size(); ++i) {
      std::cout << " t" << m_tasks[i] << " [" << est(i) << ".." << lct(i)
                << "]";
    }
    std::cout << ")\n";
  }
  size_t debug_limit{0};
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

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {

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
                << lst(i) << "]\n";
    }
    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
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

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {
      std::cout << "add t" << m_tasks[a] << ": [" << est(a) << ".."
                << minduration(a) << ".." << lct(a) << "] : " << TT.getBound()
                << std::endl;
      if (debug_flag > 1)
        std::cout << TT << std::endl;
    }
#endif

    if (TT.getBound() > lct(a)) {
      boundExplanation(lct_order.begin(), ai + 1);
      throw Failure({this, static_cast<hint>(explanations.size() - 1)});
    }
  }

  if (TT.getBound() > m_schedule.lower(HORIZON)) {
    boundExplanation(lct_order.begin(), lct_order.end());
    m_schedule.set({LOWERBOUND(HORIZON), -TT.getBound()},
                   {this, static_cast<hint>(explanations.size() - 1)});
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "]\n";
    }
    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
    }
    std::cout << " forward edges:\n";
  }
#endif

  for (auto ai{lct_order.rbegin()}; ai != (lct_order.rend() - 1); ++ai) {
    auto a{*ai};
    auto deadline_omega{lct(*(ai + 1))};
    TT.paint_gray(theta_rank[a], a);

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {
      std::cout << "rm t" << m_tasks[a] << ": [" << est(a) << ".."
                << minduration(a) << ".." << lct(a) << "] : " << TT.getBound()
                << " (" << ect(a) << ")" << std::endl;
      if (debug_flag > 1)
        std::cout << TT << std::endl;
    }
#endif

    auto ect_{TT.grayBound()};
    if (TT.getBound() > deadline_omega) {
      boundExplanation(ai + 1, lct_order.rend());
      throw Failure({this, static_cast<hint>(explanations.size() - 1)});
    }

#ifdef DEBUG_CONSTRAINT
    debug_limit = m_tasks.size();
#endif

    while (ect_ > deadline_omega) {

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << " tasks";
        for (auto j{lct_order.begin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t" << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << " can't end before " << ect_ << ", but only t"
                  << m_tasks[r] << " can go past that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r) << ".." << lct(r)
                  << "] >= " << ect_ << std::endl;
      }
#endif

      boundExplanation(ai + 1, lct_order.rend());

      explanations.emplace_back(cur_explanation);
      cur_explanation.clear();

      auto l{m_tasks[r]};
      for (auto j{lct_order.begin()}; *j != a; ++j) {

#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 0) {
          std::cout << "edge prec "
                    << m_schedule.prettyLiteral(EDGE(disjunct[r][*j]))
                    << std::endl;
        }
#endif

        //        auto t{m_tasks[*j]};
        m_schedule.set(disjunct[r][*j],
                       {this, static_cast<hint>(explanations.size() - 1)});
      }

      //      exit(1);

      if (ect(r) < ect_) {

        BoundConstraint<T> bc{LOWERBOUND(END(l)), -ect_};
#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 0) {
          std::cout << "lower bound " << bc << std::endl;
        }
#endif

        m_schedule.set(bc, {this, static_cast<hint>(explanations.size() - 1)});
      }

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << " because";
        for (auto l : explanations.back()) {
          std::cout << " ";
          m_schedule.prettyLiteral(l);
        }
        std::cout << std::endl;
      }
#endif

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 1) {
        std::cout << TT << std::endl;
        if (--debug_limit <= 0)
          exit(1);
        std::cout << debug_limit << std::endl;
      }
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

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "]\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
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

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {
      std::cout << "add t" << m_tasks[a] << ": [" << (horizon - lct(a)) << ".."
                << minduration(a) << ".." << (horizon - est(a))
                << "] : " << TT.getBound() << " (" << (horizon - lst(a)) << ")"
                << std::endl;
      if (debug_flag > 1)
        std::cout << TT << std::endl;
    }
#endif
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << m_schedule << "\nfST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "]\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
    }

    std::cout << " backward edges:\n";
  }
#endif

  for (auto ai{est_order.begin()}; ai != est_order.end() - 1; ++ai) {
    auto a{*ai};

    auto deadline_omega{horizon - est(*(ai + 1))};

    TT.paint_gray(theta_rank[a], a);

#ifdef DEBUG_CONSTRAINT
    if (debug_flag > 0) {
      std::cout << "rm t" << m_tasks[a] << ": [" << (horizon - lct(a)) << ".."
                << minduration(a) << ".." << (horizon - est(a))
                << "] : " << TT.getBound() << " (" << (horizon - lst(a)) << ")"
                << std::endl;
      if (debug_flag > 1)
        std::cout << TT << std::endl;
    }
#endif

    auto ect_{TT.grayBound()};
    if (deadline_omega < TT.getBound()) {
      boundExplanation(ai + 1, est_order.end());
      throw Failure({this, static_cast<hint>(explanations.size() - 1)});
    }

#ifdef DEBUG_CONSTRAINT
    debug_limit = m_tasks.size();
#endif

    while (ect_ > deadline_omega) {

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << "tasks";
        for (auto j{est_order.rbegin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t" << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << " can't start after " << (horizon - ect_) << ", but only t"
                  << m_tasks[r] << " can start at that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r) << ".." << lct(r)
                  << "] <= " << (horizon - ect_) << std::endl;
      }
#endif

      auto f{m_tasks[r]};
      boundExplanation(ai + 1, est_order.end());

      explanations.emplace_back(cur_explanation);
      cur_explanation.clear();

      for (auto j{est_order.rbegin()}; *j != a; ++j) {

#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 0) {
          std::cout << "edge prec "
                    << m_schedule.prettyLiteral(EDGE(disjunct[r][*j]))
                    << std::endl;
        }
#endif

        ++m_schedule.num_cons_prunings;
        m_schedule.set(disjunct[*j][r],
                       {this, static_cast<hint>(explanations.size() - 1)});
      }

      //      exit(1);

      if (lst(r) > (horizon - ect_)) {

        BoundConstraint<T> bc{UPPERBOUND(START(f)), (horizon - ect_)};
#ifdef DEBUG_CONSTRAINT
        if (debug_flag > 0) {
          std::cout << "upper bound " << bc << std::endl;
        }
#endif

        m_schedule.set({UPPERBOUND(START(f)), (horizon - ect_)},
                       {this, static_cast<hint>(explanations.size() - 1)});
      }

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 0) {
        std::cout << " because";
        for (auto l : explanations.back()) {
          std::cout << " ";
          m_schedule.prettyLiteral(l);
        }
        std::cout << std::endl;
      }
#endif

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DEBUG_CONSTRAINT
      if (debug_flag > 1) {
        std::cout << TT << std::endl;
        if (--debug_limit <= 0)
          exit(1);
        std::cout << debug_limit << std::endl;
      }
#endif
    }
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i) << ".."
                << lst(i) << "]\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
    }
  }
#endif
}

template <typename T> int DisjunctiveEdgeFinding<T>::getType() const {
  return EDGEFINDINGEXPL;
}

template <typename T>
void DisjunctiveEdgeFinding<T>::xplain(const lit, const hint h,
                                       std::vector<lit> &Cl) {
  for (auto r : explanations[h]) {
    Cl.push_back(r);
  }
}

template <typename T>
std::ostream &DisjunctiveEdgeFinding<T>::display(std::ostream &os) const {
  os << "Disjunctive Edge-Finding";

#ifdef DEBUG_CONSTRAINT
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
