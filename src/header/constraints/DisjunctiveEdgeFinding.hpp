/************************************************
 * Tempo DisjunctiveEdgeFinding.hpp
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

#ifndef TEMPO_DISJUNCTIVEEDGEFINDING_HPP
#define TEMPO_DISJUNCTIVEEDGEFINDING_HPP

#include <cassert>
#include <map>
#include <vector>
#include <ranges>
#include <Iterators.hpp>

#include "Explanation.hpp"
#include "Global.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"
#include "util/ThetaTree.hpp"
#include "util/Matrix.hpp"
#include "util/traits.hpp"
#include "Model.hpp"


namespace tempo {

template<typename T>
class Solver;

template <typename T> class DisjunctiveEdgeFinding : public Constraint<T> {
private:
  Solver<T> &m_solver;
  Interval<T> schedule;
  std::vector<Interval<T>> the_tasks;
  Matrix<Literal<T>> disjunct;
    
    std::vector<T> est_buffer;
    std::vector<T> lct_buffer;
    bool backward_flag{false};

  // helpers
  std::vector<unsigned> est_order;
  std::vector<unsigned> lct_order;
  std::vector<unsigned> theta_rank;
  ThetaTree TT;

//  std::vector<std::pair<Literal<T>, hint>> pruning;
  std::vector<std::vector<Interval<T> *>> the_explanation_tasks;
  std::vector<T> explanation_lb;
  std::vector<T> explanation_ub;

  Reversible<size_t> num_explanations;

  std::vector<unsigned> pruned_tasks;
  std::vector<std::vector<unsigned>::reverse_iterator> omegas;
  std::vector<std::vector<unsigned>::iterator> bomegas;
  std::vector<T> relevant_starts;
  std::vector<T> bound_omegas;

  int level;
    
#ifdef DBG_EDGEFINDING
    bool checkpruning(const Literal<T> prec, const hint h);
    bool checkoverload(const T lb, const hint h);
#endif
    
  template <typename Iter>
  hint boundExplanation(Iter b, Iter e, const T lb);

  void propagateForward();
    
    void getDomains();
    void negateDomains();
    void applyPruning();

public:
  template <concepts::typed_range<Interval<T>> Tasks, typename ItVar> requires(std::ranges::sized_range<Tasks>)
  DisjunctiveEdgeFinding(Solver<T> &solver, Interval<T> &sched, Tasks &&tasks, ItVar beg_var);
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
//  int getType() const override;

    std::string asciiArt(const int i) const;
  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;

  void printLBExplanation(const hint ph);
  void printUBExplanation(const hint ph);
  void printTrivialExplanation(const Literal<T> l);

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
  ss << "t" << the_tasks[i].id() << ": [" << est_buffer[i] << ".." << lct_buffer[i] << "] ("
     << minduration(i) << ")";
  return ss.str();
}

#ifdef DBG_EDGEFINDING
template <typename T>
bool DisjunctiveEdgeFinding<T>::checkpruning(const Literal<T> l, const hint h) {
    
    T total_duration{0};
    
    for (auto ti : the_explanation_tasks[h]) {
        total_duration += ti->minDuration(m_solver);
    }
    
    return total_duration > (explanation_ub[h] - explanation_lb[h]);
    
//    if (l.isNumeric()) {
//      // explain the bound change
//      if (l.sign() == bound::lower) {
//          return explanation_lb[h] + total_duration >= l.value();
//      } else {
//          return explanation_ub[h] - total_duration <= -l.value();
//      }
//    } else {
//      // explain the edges
//      auto lc{m_solver.boolean.getEdge(l)};
//      if (t->getStart() == lc.from or t->getEnd() == lc.from) {
//        Cl.push_back(t->start.after(explanation_lb[h]));
//      } else {
//        Cl.push_back(t->end.before(explanation_ub[h]));
//      }
//    }
}

template <typename T>
bool DisjunctiveEdgeFinding<T>::checkoverload(const T lb, const hint h) {
    auto ect{explanation_lb[h]};
    
    for (auto ti : the_explanation_tasks[h]) {
        ect += ti->minDuration(m_solver);
    }
    
    return ect >= lb;
}
#endif

// collect all the tasks in [b,e) whose lower bound is larger than or equal to
// lb
template <typename T>
template <typename Iter>
hint DisjunctiveEdgeFinding<T>::boundExplanation(Iter b, Iter e,
                                                      const T lb) {
  auto e_idx{num_explanations};
  if (the_explanation_tasks.size() <= e_idx) {
    the_explanation_tasks.resize(e_idx + 1);
    explanation_lb.resize(e_idx + 1);
    explanation_ub.resize(e_idx + 1);
  } else {
    the_explanation_tasks[e_idx].clear();
  }
    
    if(backward_flag)
        explanation_ub[e_idx] = -lb;
    else
        explanation_lb[e_idx] = lb;

  T ub{lct_buffer[*b]};
  for (auto x{b}; x != e; ++x) {
    if (est_buffer[*x] >= lb) {

#ifdef DBG_EXPLEF
      std::cout << " expl: " the_tasks[*x].start.after(lb);
<<)
<< " and "
the_tasks[*x].end.before(the_tasks[*x].getLatestEnd(m_solver));
<< std::endl;
#endif

//      explanation_tasks[e_idx].push_back(i);
the_explanation_tasks[e_idx].push_back(&the_tasks[*x]);

ub = std::max(ub, lct_buffer[*x]);
    }
#ifdef DBG_EXPLEF
    else {
      std::cout << " skip " << prettyTask(*x) << " b/c of the lower bound ("
                << lb << ")" << std::endl;
    }
#endif
  }

  assert(ub >= lct_buffer[*b]);
    
    if(backward_flag)
        explanation_lb[e_idx] = -ub;
    else
        explanation_ub[e_idx] = ub;

  ++num_explanations;

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
  return the_tasks[i].minDuration(m_solver);
}

template <typename T>
T DisjunctiveEdgeFinding<T>::maxduration(const unsigned i) const {
  return the_tasks[i].maxDuration(m_solver);
}

template <typename T>
template <concepts::typed_range<Interval<T>> Tasks, typename ItVar> requires(std::ranges::sized_range<Tasks>)
DisjunctiveEdgeFinding<T>::DisjunctiveEdgeFinding(Solver<T> &solver, Interval<T> &sched, Tasks &&tasks,
                                                  ItVar beg_var) : m_solver(solver), schedule(sched),
                                                                   the_tasks(std::forward<Tasks>(tasks).begin(),
                                                                             std::forward<Tasks>(tasks).end()),
                                                                   disjunct(the_tasks.size(), the_tasks.size()),
                                                                   est_buffer(the_tasks.size()),
                                                                   lct_buffer(the_tasks.size()), TT(the_tasks.size()),
                                                                   num_explanations(0, &(m_solver.getEnv())) {
  using iterators::const_enumerate;
  using std::views::drop;
  Constraint<T>::priority = Priority::Medium;

  for (unsigned i = 0; i < the_tasks.size(); ++i) {
    lct_order.push_back(i);
    est_order.push_back(i);
  }

  auto ep{beg_var};
  for(auto [i, ip]: const_enumerate(the_tasks)) {
    for (auto [j, jp]: const_enumerate(the_tasks | drop(i + 1), i + 1)) {
      auto x{*ep};
      disjunct(i, j) = m_solver.boolean.getLiteral(true, x);
      disjunct(j, i) = m_solver.boolean.getLiteral(false, x);
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
bool DisjunctiveEdgeFinding<T>::notify(const Literal<T>, const int) {
  return true;
}

template <typename T> void DisjunctiveEdgeFinding<T>::getDomains() {
    for(unsigned i{0}; i<the_tasks.size(); ++i) {
        est_buffer[i] = est(i);
        lct_buffer[i] = lct(i);
    }
    backward_flag = false;
}

template <typename T> void DisjunctiveEdgeFinding<T>::negateDomains() {
    for(unsigned i{0}; i<the_tasks.size(); ++i) {
        est_buffer[i] = -lct(i);
        lct_buffer[i] = -est(i);
    }
    backward_flag = true;
}

//template <typename T> void DisjunctiveEdgeFinding<T>::applyPruning() {
//}

//template <typename T> void DisjunctiveEdgeFinding<T>::negateExplanations(size_t &ptr) {
//    //
//    auto end_exp{static_cast<size_t>(num_explanations)};
//    while(ptr < end_exp) {
//        explanation_lb[ptr] = -explanation_lb[ptr];
//    }
//    
//    for(size_t i{0}; i<pruning.size(); ++i) {
//        
//    }
//}
              
//              template <typename T> void DisjunctiveEdgeFinding<T>::applyPruning() {
//                  for(auto p : pruning) {
//                      m_solver.set(p.first, {this, p.second});
//                  }
//              }

template <typename T> void DisjunctiveEdgeFinding<T>::propagate() {
    
    getDomains();
    
#ifdef DBG_EDGEFINDING
  if (DBG_EDGEFINDING) {
        std::cout << "\n\npropagate edge-finding  (";
        std::cout << prettyTask(0);
        for (size_t i{1}; i < the_tasks.size(); ++i) {
            std::cout << " " << prettyTask(i);
        }
        std::cout << ") (i=" << m_solver.num_cons_propagations << ")\n";
      }
    #endif
    
    
  propagateForward();
    applyPruning();
    
    negateDomains();
    propagateForward();
    applyPruning();

}

template <typename T> void DisjunctiveEdgeFinding<T>::propagateForward() {
    
    //
    //    std::cout << std::endl << "propagation pass (";
    //    std::cout << prettyTask(0);
    //    for (size_t i{1}; i < the_tasks.size(); ++i) {
    //        std::cout << " " << prettyTask(i);
    //    }
    //    std::cout << ") (i=" << m_solver.num_cons_propagations << ")\n";
    //  }
    //#endif
    
    std::sort(est_order.begin(), est_order.end(),
              [&](const unsigned a, const unsigned b) -> bool {
        return est_buffer[a] < est_buffer[b];
    });
    
    std::sort(lct_order.begin(), lct_order.end(),
              [&](const unsigned a, const unsigned b) -> bool {
        return lct_buffer[a] < lct_buffer[b]; // or (lct_buffer[a] == lct_buffer[b] and est_buffer[a] < est_buffer[b]);
    });
    
#ifdef DBG_EDGEFINDING
    if (DBG_EDGEFINDING) {
        std::cout << "\npropagation pass\n";
        for (auto j : lct_order) {
            std::cout << "task t" << the_tasks[j].id() << ": " << asciiArt(j) << std::endl;
        }
    }
#endif
    
    TT.clear();
    
    for (unsigned i{0}; i < est_order.size(); ++i) {
        theta_rank[est_order[i]] = i;
    }
    
    // insert the tasks one by one in the theta tree by non-decreasing lct order
    for (auto ai{lct_order.begin()}; ai != lct_order.end(); ++ai) {
        auto a{*ai};
        TT.insert(theta_rank[a], est_buffer[a], minduration(a));
        
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
            std::cout << "insert " << prettyTask(a) << " bound=" << TT.getBound()
             << ", relevant est=" << TT.getEst() << std::endl;
        }
#endif
        
        // failure because of overload
        if (TT.getBound() > lct_buffer[a]) {
            auto h{boundExplanation(lct_order.begin(), ai + 1, TT.getEst())};
            
#ifdef DBG_EDGEFINDING
            if (DBG_EDGEFINDING) {
                std::cout << " failure ect_omega = " << TT.getBound() << " > lct_omega = " << lct_buffer[a]
                << std::endl;
            }
            
            if(not checkoverload(TT.getBound(),h)) {
                std::cout << "bug overload fail\n";
                exit(1);
            }
#endif
            
            throw Failure<T>({this, h});
        }
    }
    
//    std::cout << "schedule.getEarliestEnd(m_solver) = " << schedule.getEarliestEnd(m_solver) << std::endl;
    
    // global upper bound from the overload formula
    if (TT.getBound() > schedule.getEarliestEnd(m_solver)) {
        
        auto h{boundExplanation(lct_order.begin(), lct_order.end(), TT.getEst())};
//        pruning.emplace_back(schedule.end.after(-TT.getBound()), h);
        
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
            std::cout << " makespan pruning: " << m_solver.pretty(schedule.end.after(TT.getBound())) << " b/c" ;
            for(auto t : the_explanation_tasks[h]) {
                std::cout << " t" << t->id();
            }
            std::cout << " in [" << explanation_lb[h] << ".." << explanation_ub[h] << "]\n"; // << std::endl;
        }
        
        if(not checkoverload(TT.getBound(),h)) {
            std::cout << "bug overload pruning\n";
            exit(1);
        }
#endif
        
        m_solver.set(schedule.end.after(TT.getBound()), {this, h});
        
//        std::cout << "schedule.getEarliestEnd(m_solver) = " << schedule.getEarliestEnd(m_solver) << std::endl;
        
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
        auto deadline_omega{lct_buffer[*(ai + 1)]};
        TT.paint_gray(theta_rank[a], a);
        
        // bound including the processing of the gray task
        auto ect_{TT.grayBound()};
        assert(TT.getBound() <= deadline_omega);
        
#ifdef DBG_EDGEFINDING
        if (DBG_EDGEFINDING) {
            std::cout << "gray " << prettyTask(a) << " lct=" << deadline_omega
            << ", bounds=[" << TT.grayEst() << ".." << ect_ << "]\n";// << std::endl;
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
}


template <typename T> void DisjunctiveEdgeFinding<T>::applyPruning() {

  hint ph{Constant::NoHint};

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
        ph = Constant::NoHint;
        for (auto j{ai}; j != lct_order.rend(); ++j) {
            
            auto prec{(backward_flag ? disjunct(*j, r) : disjunct(r, *j))};
            
#ifdef DBG_EDGEFINDING
            if (DBG_EDGEFINDING) {
                std::cout << "add precedences " << m_solver.pretty(prec) << "?" << std::endl;
            }
#endif
            
            if (not m_solver.boolean.satisfied(prec)) {
                
                // ej < si (ub(si) & lb(ej))
//                if (m_solver.boolean.falsified(prec)) {
//                    // the precedence is trivially implied because r >> *j
//                    // not clear if we should just let the edge constraints handle that
//                    
//                    // this is the edge that is satisfied and implies disjunct[r][*j]
//                    auto eij{~m_solver.boolean.getEdge(~prec)};
//                    auto idx_lb{m_solver.numeric.lastLitIndex(bound::lower, eij.from)};
//                    auto idx_ub{m_solver.numeric.lastLitIndex(bound::upper, eij.to)};
//                    
//                    
//#ifdef DBG_EDGEFINDING
//                    if (DBG_EDGEFINDING) {
//                        std::cout << " trivial precedence: " << m_solver.pretty(prec) << std::endl;
//                    }
//                    
////                    if(not checktrivialpruning(prec)) {
////                        std::cout << "bug trivial pruning\n";
////                        exit(1);
////                    }
//#endif
//                    
//                    // ei < sj (ub(ei) & lb(sj))
//                    //          pruning.emplace_back(disjunct[r][*j], -1 - static_cast<hint>(std::min(idx_lb, idx_ub)));
//                    m_solver.set(prec, {this, -1 - static_cast<hint>(std::min(idx_lb, idx_ub))});
//                    
//                } else {
                
                // the precedences where the precedence is already falsified will be handled by disjunctions
                if (not m_solver.boolean.falsified(prec)) {
                  if (ph == Constant::NoHint) {
                    ph = boundExplanation(ai, lct_order.rend(), s);
                    the_explanation_tasks[ph].push_back(&tl);
                  }

#ifdef DBG_EDGEFINDING
                    if (DBG_EDGEFINDING) {
                        std::cout << " edge: " << m_solver.pretty(prec) << " b/c" ;
                        for(auto t : the_explanation_tasks[ph]) if(t->id() != tl.id()) {
                            std::cout << " t" << t->id();
                        }
                        std::cout << " in [" << explanation_lb[ph] << ".." << explanation_ub[ph] << "] and t" << tl.id();
                        if(backward_flag)
                            std::cout << " end before " << explanation_ub[ph] << "\n";
                        else
                            std::cout << " start after " << explanation_lb[ph] << "\n";
                    }
                    
                    if(not checkpruning(prec,ph)) {
                        std::cout << "bug prec\n";
                        exit(1);
                    }
#endif
                    
                    //          pruning.emplace_back(disjunct[r][*j], ph);
                    m_solver.set(prec, {this, ph});
                    
                    
                }
            }
        }
        
        // not est because 1/ ect_omega includes task tl, and 2/ it is still correct in the preemptive case
        bool relevant{(backward_flag ? lst(r) > -ect_omega : ect(r) < ect_omega)};
        if (relevant) {
            
            Literal<T> bc{(backward_flag ? tl.start.before(-ect_omega) : tl.end.after(ect_omega))};

            if (ph == Constant::NoHint) {
              ph = boundExplanation(ai, lct_order.rend(), s);
              the_explanation_tasks[ph].push_back(&tl);
            }

#ifdef DBG_EDGEFINDING
            if (DBG_EDGEFINDING) {
                std::cout << " adjustment: " << m_solver.pretty(bc) << " b/c" ;
                for(auto t : the_explanation_tasks[ph]) if(t->id() != tl.id()) {
                    std::cout << " t" << t->id();
                }
                std::cout << " in [" << explanation_lb[ph] << ".." << explanation_ub[ph] << "] and t" << tl.id() ;
                if(backward_flag)
                    std::cout << " end before " << explanation_ub[ph] << "\n";
                else
                    std::cout << " start after " << explanation_lb[ph] << "\n";
            }
#endif
            
            //        pruning.emplace_back(bc, ph);
            m_solver.set(bc, {this, ph});
            
            
        }
        
#ifdef DBG_EDGEFINDING
        else if (DBG_EDGEFINDING) {
            std::cout << "already satisfied" << std::endl;
        }
#endif
    }
}



template <typename T>
void DisjunctiveEdgeFinding<T>::xplain(const Literal<T> l, const hint h,
                                       std::vector<Literal<T>> &Cl) {

  //  if (static_cast<size_t>(h) >= explanation_tasks.size()) {
  //    std::cout << h << " / " << explanation_tasks.size() << std::endl;
  //    exit(1);
  //  }

  if (l == Contradiction<T>) {
#ifdef DBG_EXPLEF
    std::cout << "explain failure from edge-finding: overload on interval ["
              << explanation_lb[h] << ".." << explanation_ub[h] << "]\n";
    T duration{0};
#endif

    // failure case, everything is in "tasks"
    assert(static_cast<std::size_t>(h) < the_explanation_tasks.size());
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
std::string DisjunctiveEdgeFinding<T>::asciiArt(const int i) const {
  std::stringstream ss;
  ss << std::setw(3)
     << std::left << minduration(i) << " " << std::right;
  for (auto k{0}; k < est(i); ++k) {
    ss << " ";
  }
  ss << "[";
  for (auto k{est(i) + 1}; k < ect(i); ++k) {
    ss << "=";
  }
    if(lct(i) == Constant::Infinity<T>) {
        ss << "=... " << est(i) << "..\\inf";
    } else {
        for (auto k{ect(i)}; k < lct(i) - 1; ++k) {
            ss << ".";
        }
        ss << "] " << est(i) << ".." << lct(i);
    }
  return ss.str();
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
