#ifndef TEMPO_DISJUNCTIVEEDGEFINDING_HPP
#define TEMPO_DISJUNCTIVEEDGEFINDING_HPP

#include <cassert>
#include <map>
#include <vector>

#include "Global.hpp"
#include "Constraint.hpp"
#include "Scheduler.hpp"
#include "Explanation.hpp"
#include "utilSparseSet.hpp"
#include "util/ThetaTree.hpp"

namespace tempo {

template <class T>
class DisjunctiveEdgeFinding : public Constraint, public Explainer {
private:
  Scheduler<T> &m_schedule;
  std::vector<task> m_tasks;

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
                         const ItTask end_task,
                         const ItVar beg_var,
                         const ItVar end_var
                         );
  virtual ~DisjunctiveEdgeFinding();

  // helpers
  T est(const unsigned i) const;
  T lst(const unsigned i) const;
  T ect(const unsigned i) const;
  T lct(const unsigned i) const;
  T minduration(const unsigned i) const;
  T maxduration(const unsigned i) const;

  bool notify(const int var, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
  void xplain(const lit l, const hint h, CutBuilder *CB) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  C getCapacity() const noexcept { return 1; }

  auto getTasks() const noexcept -> const std::vector<task> & {
    return m_tasks;
  }

  auto getTaskConsumptions() const noexcept -> const std::vector<C> & {
    std::vector cons(m_tasks.size(), 1);
    return cons;
  }
    
#ifdef DEBUG_CONSTRAINT
    debug_flag{false};
#endif
};

template <class T, class C>
template <typename Iter>
void DisjunctiveEdgeFinding<T, C>::boundExplanation(Iter b, Iter e) {

  assert(cur_explanation.empty());

  for (auto x{b}; x != e; ++x) {

    auto i{m_tasks[*x]};

    distance.pushLiteralsOn(START(i), ORIGIN, cur_explanation);
    distance.pushLiteralsOn(ORIGIN, START(i), cur_explanation);
    distance.pushLiteralsOn(END(i), START(i), cur_explanation);
  }

  explanations.push_back(cur_explanation);
//  ++m_schedule.num_cons_prunings;
  cur_explanation.clear();
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::est(const unsigned i) const {
  return m_scheduler.lower(START(m_tasks[i]));
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::lst(const unsigned i) const {
  return m_scheduler.upper(START(m_tasks[i]));
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::ect(const unsigned i) const {
  return m_scheduler.lower(END(m_tasks[i]));
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::lct(const unsigned i) const {
  return m_scheduler.upper(END(m_tasks[i]));
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::minduration(const unsigned i) const {
  return m_scheduler.minDuration(m_tasks[i]);
}

template <class T, class C>
T DisjunctiveEdgeFinding<T, C>::maxduration(const unsigned i) const {
    return m_scheduler.maxDuration(m_tasks[i]);
}

template <class T, class C>
template <typename ItTask, typename ItVar>
DisjunctiveEdgeFinding<T, C>::DisjunctiveEdgeFinding(Scheduler<T> &scheduler,
                                                     const ItTask beg_task,
                                                     const ItTask end_task,
                                                     const ItVar beg_var,
                                                     const ItVar end_var)
: m_schedule(scheduler), TT(std::distance(beg_task, end_task)) {
    
    // get all tasks with non-zero duration
    for (auto j{beg_task}; j != end_task; ++j) {
        task t{*j};
        if (maxduration(t) > 0) {
            m_tasks.push_back(t);
        }
    }
    
    disjunct.resize(m_tasks.size());
    
    for (unsigned i = 0; i < m_tasks.size(); ++i) {
        lct_order.push_back(i);
        est_order.push_back(i);
        disjunct[i].resize(m_tasks.size());
    }
    
    for(auto v{beg_var}; v!=end_var; ++v) {
        auto e{m_schedule.getEdge(*v)};
        auto i{TASK(e.from)};
        auto j{TASK(e.to)};
        disjunct[i][j] = POS(var);
        disjunct[j][i] = NEG(var);
    }
    
    theta_rank.resize(m_tasks.size(), 0);
}

template <class T, class C>
DisjunctiveEdgeFinding<T, C>::~DisjunctiveEdgeFinding() {}

template <class T, class C>
void DisjunctiveEdgeFinding<T, C>::post(const int idx) {

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  cons_id = idx;

    for (unsigned i{0}; i < m_tasks.size(); ++i) {
        m_schedule.wake_me_on_event(LOWERBOUND(START(m_tasks[i])), cons_id);
        m_schedule.wake_me_on_event(UPPERBOUND(END(m_tasks[i])), cons_id);
    }

}

template <class T, class C>
bool DisjunctiveEdgeFinding<T, C>::notify(const int var, const int rank) {
    return true;
}

template <class T, class C> void DisjunctiveEdgeFinding<T, C>::propagate() {

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {
    std::cout << std::endl << "propagate edge-finding(";
    for (auto t : m_tasks) {
      std::cout << " t" << t << " [" << est(t) << ".." << lct(t) << "]";
    }
    std::cout << ")\n";
  }
  size_t size_before{0};
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

    auto m_schedule.upper(HORIZON);

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {

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
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
                << ".." << lst(i) << "]\n";
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
    if (this->ResourceConstraint<T>::debug_flag) {
      std::cout << "add " << m_schedule.taskLabel(m_tasks[a]) << ": [" << est(a)
                << ".." << minduration(a) << ".." << lct(a)
                << "] : " << TT.getBound() << std::endl;
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
      m_schedule.set({LOWERBOUND(HORIZON), -TT.getBound()}, {this, static_cast<hint>(explanations.size() - 1)});
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
                << ".." << lst(i) << "]\n";
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
    if (this->ResourceConstraint<T>::debug_flag) {
      std::cout << "rm t" << m_tasks[a] << ": [" << est(a)
                << ".." << minduration(a) << ".." << lct(a)
                << "] : " << TT.getBound() << " (" << ect(a) << ")"
                << std::endl;
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
      if (debug_flag) {
        std::cout << " tasks";
        for (auto j{lct_order.begin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t"
                  << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DEBUG_CONSTRAINT
      if (this->ResourceConstraint<T>::debug_flag) {
        std::cout << " can't end before " << ect_ << ", but only t"
                  << m_tasks[r]
                  << " can go past that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r)
                  << ".." << lct(r) << "] >= " << ect_ << std::endl;
      }
#endif

      boundExplanation(ai + 1, lct_order.rend());
      auto l{m_tasks[r]};
      for (auto j{lct_order.begin()}; *j != a; ++j) {
        auto t{m_tasks[*j]};
          m_schedule.set(disjunct[l][t],{this, static_cast<hint>(explanations.size() - 1)});
      }

      if (ect(r) < ect_) {
          m_schedule.set({LOWERBOUND(END(l)), -ect_}, {this, static_cast<hint>(explanations.size() - 1)});
      }

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DEBUG_CONSTRAINT
      if (debug_flag) {
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
  if (debug_flag) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
                << ".." << lst(i) << "]\n";
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
    if (debug_flag) {
      std::cout << "add t" << m_tasks[a] << ": ["
                << (horizon - lct(a)) << ".." << minduration(a) << ".."
                << (horizon - est(a)) << "] : " << TT.getBound() << " ("
                << (horizon - lst(a)) << ")" << std::endl;
      std::cout << TT << std::endl;
    }
#endif
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {
    std::cout << std::endl << m_schedule << "\nfST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
                << ".." << lst(i) << "]\n";
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
    if (debug_flag) {
      std::cout << "rm t" << m_tasks[a] << ": ["
                << (horizon - lct(a)) << ".." << minduration(a) << ".."
                << (horizon - est(a)) << "] : " << TT.getBound() << " ("
                << (horizon - lst(a)) << ")" << std::endl;
      std::cout << TT << std::endl;
    }
#endif

    auto ect_{TT.grayBound()};
    if(deadline_omega < TT.getBound()) {
      boundExplanation(ai + 1, est_order.end());
      throw Failure({this, static_cast<hint>(explanations.size() - 1)});
    }


#ifdef DEBUG_CONSTRAINT
    debug_limit = m_tasks.size();
#endif

    while (ect_ > deadline_omega) {

#ifdef DEBUG_CONSTRAINT
      if (debug_flag) {
        std::cout << "tasks";
        for (auto j{est_order.rbegin()}; *j != a; ++j)
          std::cout << " t" << m_tasks[*j];
        std::cout << " + t"
                  << m_tasks[TT.getResponsible()];
      }
#endif

      auto r{TT.getResponsible()};

#ifdef DEBUG_CONSTRAINT
      if (debug_flag) {
        std::cout << " can't start after " << (horizon - ect_) << ", but only t"
                  << m_tasks[r]
                  << " can start at that time\n ==> t"
                  << m_tasks[r] << " in [" << est(r)
                  << ".." << lct(r) << "] <= " << (horizon - ect_) << std::endl;
      }
#endif

      auto f{m_tasks[r]};
      boundExplanation(ai + 1, est_order.end());
      for (auto j{est_order.rbegin()}; *j != a; ++j) {
        auto t{m_tasks[*j]};

        ++m_schedule.num_cons_prunings;
          m_schedule.set(disjunct[t][f], {this, static_cast<hint>(explanations.size() - 1)});
      }

      if (lst(r) > (horizon - ect_)) {

        explanations.push_back(cur_explanation);
//        ++m_schedule.num_cons_prunings;
        cur_explanation.clear();

          m_schedule.set({UPPERBOUND(START(f)), (horizon - ect_)}, {this, static_cast<hint>(explanations.size() - 1)});
      }

      TT.remove(theta_rank[r]);

      ect_ = TT.grayBound();

#ifdef DEBUG_CONSTRAINT
      if (debug_flag) {
        std::cout << TT << std::endl;
        if (--debug_limit <= 0)
          exit(1);
        std::cout << debug_limit << std::endl;
      }
#endif
    }
  }

#ifdef DEBUG_CONSTRAINT
  if (debug_flag) {
    std::cout << std::endl << "ST\n";
    for (auto i : est_order) {
      std::cout << prettyEvent(START(m_tasks[i])) << ": [" << est(i)
                << ".." << lst(i) << "]\n";
    }

    std::cout << "\nCT\n";
    for (auto i : lct_order) {
      std::cout << prettyEvent(END(m_tasks[i])) << ": [" << ect(i) << ".."
                << lct(i) << "]\n";
    }
  }
#endif

}

template <class T, class C> int DisjunctiveEdgeFinding<T, C>::getType() const {
  return EDGEFINDINGEXPL;
}

template <class T, class C>
void DisjunctiveEdgeFinding<T, C>::xplain(const lit l, const hint h, std::vector<lit>& Cl) {
  for (auto r : explanations[h]) {
    Cl.push_back(r);
  }
}

template <class T, class C>
std::ostream &DisjunctiveEdgeFinding<T, C>::display(std::ostream &os) const {
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

template <class T, class C>
std::ostream &DisjunctiveEdgeFinding<T, C>::print_reason(std::ostream &os,
                                                         const hint h) const {
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

} // namespace schedcl

#endif
