/************************************************
 * Tempo CumulativeCheck.hpp
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

#ifndef TEMPO_CUMULATIVE_HPP
#define TEMPO_CUMULATIVE_HPP

#include <cassert>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "constraints/Constraint.hpp"
#include "util/DisjointSet.hpp"
#include "util/SparseSet.hpp"
#include "util/LexBFS.hpp"
#include "Model.hpp"

//#define CHECK_CLIQUE

namespace tempo {

template<typename T>
class Solver;

template <typename T> class CumulativeCheck : public Constraint<T> {
private:
  Solver<T> &m_solver;
  NumericVar<T> capacity;
  std::vector<Interval<T>> the_tasks;
  std::vector<NumericVar<T>> demand;

  // tasks that are in relevant with at least another task
  SparseSet<int, Reversible<size_t>> relevant;

  // j in relevant[i] <=> arc (i,j) <=> task[i] and task[j] are in relevant
  std::vector<SparseSet<int, Reversible<size_t>>> parallel;
  // precedence[i][j] <=> (e_i <= s_j or e_i > s_j)
  std::vector<std::vector<Literal<T>>> precedence;
  std::vector<int> scopex;
  std::vector<int> scopey;

  LexBFS simplicial;

  //    std::vector<unsigned> sorted_tasks;

  std::vector<int> fail_xpl;
  SparseSet<> clique;
    
    bool pruning_flag{false};

public:
  template <typename ItTask, typename ItNVar, typename ItBVar>
  CumulativeCheck(Solver<T> &solver, const NumericVar<T> c,
                  const ItTask beg_task, const ItTask end_task,
                  const ItNVar beg_dem, const ItBVar beg_disj);
  virtual ~CumulativeCheck();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;
//  int getType() const override;

  bool start_before_end(const int i, const int j) const;
  bool end_before_start(const int i, const int j) const;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

  std::string prettyTask(const int i) const;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename ItTask, typename ItNVar, typename ItBVar>
CumulativeCheck<T>::CumulativeCheck(Solver<T> &solver, const NumericVar<T> c,
                                    const ItTask beg_task,
                                    const ItTask end_task, const ItNVar beg_dem,
                                    ItBVar beg_disj)
    : m_solver(solver), capacity(c),
      relevant(std::distance(beg_task, end_task), &solver.getEnv()) {

  Constraint<T>::priority = Priority::Low;

  auto dp{beg_dem};
  for (auto jp{beg_task}; jp != end_task; ++jp) {
    the_tasks.push_back(*jp);
    demand.push_back(*dp);

    //              std::cout << *jp << " requires " << *dp <<
    //              std::endl;
    ++dp;
  }

  clique.reserve(the_tasks.size());
  precedence.resize(the_tasks.size());

  relevant.clear();

  for (size_t i{0}; i < the_tasks.size(); ++i) {

    precedence[i].resize(the_tasks.size());

    parallel.emplace_back(the_tasks.size(), &m_solver.getEnv());
    parallel.back().clear();
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

      //        std::cout << *ip << " // " << *jp << " iff (" <<
      //        precedence[i][j] << " and " << precedence[j][i] << ")\n";

      ++ep;
    }
  }
}

template <typename T> CumulativeCheck<T>::~CumulativeCheck() {}

//template <typename T> bool CumulativeCheck<T>::relevant(const int i, const int j) const {
//        return precedence[i][j].satisfied(m_solver) and precedence[j][i].satisfied(m_solver)
//}

template <typename T> bool CumulativeCheck<T>::start_before_end(const int i, const int j) const {
  return m_solver.boolean.satisfied(precedence[j][i]);
}

template <typename T>
bool CumulativeCheck<T>::end_before_start(const int i, const int j) const {
  return m_solver.boolean.falsified(precedence[i][j]);
}

template <typename T> void CumulativeCheck<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (size_t i{0}; i < the_tasks.size(); ++i) {
    for (size_t j{0}; j < the_tasks.size(); ++j)
      if (i != j) {
        m_solver.wake_me_on(precedence[i][j], this->id());
        scopex.push_back(i);
        scopey.push_back(j);
      }
  }
}

template <typename T>
bool CumulativeCheck<T>::notify(const Literal<T>
#ifdef DBG_CCHECK
                                    l
#endif
                                ,
                                const int r) {

  auto x{scopex[r]};
  auto y{scopey[r]};

#ifdef DBG_CCHECK
  if (DBG_CCHECK) {
    std::cout << "\nnotify  end(" << the_tasks[x] << ") > start("
              << the_tasks[y] << ") / " << m_solver.pretty(l) << std::endl;

    std::cout << start_before_end(y, x) << "|" << start_before_end(x, y)
              << std::endl;

    std::cout << "level=" << m_solver.level() << " relevant: " << relevant
              << std::endl;

    for (unsigned i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << ": ";
      for (unsigned j{0}; j < the_tasks.size(); ++j) {
        if (i != j and start_before_end(i, j) and start_before_end(j, i)) {
          if (parallel[i].has(j))
            std::cout << "*";
          else
            std::cout << "?";
        } else {
          if (not parallel[i].has(j))
            std::cout << ".";
          else
            std::cout << "~";
        }

        assert(parallel[i].has(j) == parallel[j].has(i));
      }
      std::cout << std::endl;
    }

    for (unsigned i{0}; i < the_tasks.size(); ++i) {
      std::cout << "t" << the_tasks[i].id() << " ["
                << the_tasks[i].getEarliestStart(m_solver) << "-"
                << the_tasks[i].getLatestStart(m_solver) << ".."
                << the_tasks[i].getEarliestEnd(m_solver) << "-"
                << the_tasks[i].getLatestEnd(m_solver)
                << "] p=" << the_tasks[i].minDuration(m_solver)
                << " d=" << demand[i].min(m_solver) << std::endl;
    }
  }
#endif

  if (start_before_end(x, y)) {
    if (not relevant.has(x)) {
      relevant.add(x);
    }
    if (not relevant.has(y)) {
      relevant.add(y);
    }

    //        assert(not parallel[x].has(y) and not parallel[y].has(x));
    if (not parallel[x].has(y)) {
      parallel[x].add(y);
      parallel[y].add(x);
    }

    return true;
  }

#ifdef DBG_CCHECK
  if (DBG_CCHECK) {
    std::cout << "no new parallelism\n";
  }
#endif

  return false;
}

template <typename T> void CumulativeCheck<T>::propagate() {

#ifdef DBG_CCHECK
  if (DBG_CCHECK) {
    std::cout << "propagate " << *this << std::endl;
    std::cout << "relevant: " << relevant << std::endl;
  }
#endif

  simplicial.explore(relevant, parallel);

#ifdef DBG_CCHECK
  if (DBG_CCHECK) {
    std::cout << "simplicial elimination:" << std::endl;
  }
#endif

  for (auto v{simplicial.ordering.rbegin()}; v != simplicial.ordering.rend();
       ++v) {

#ifdef DBG_CCHECK
    if (DBG_CCHECK) {
      std::cout << " * clique " << *v << ":";
      for (auto w : parallel[*v])
        if (simplicial.before(*v, w)) {
          std::cout << " " << w;
          for (auto u : parallel[*v])
            if (simplicial.before(*v, u) and u != w) {
              if (not parallel[u].has(w) or not parallel[w].has(u)) {
                std::cout << std::endl
                          << *v << " is not simplicial because " << u << " and "
                          << w << " are not adjacent\n";
                exit(1);
              }
            }
        }
      std::cout << std::endl;
    }
#endif

    clique.clear();
    clique.add(*v);
    T size_clique{demand[*v].min(m_solver)};
    for (auto w : parallel[*v]) {
        if (simplicial.before(*v, w)) {
            size_clique += demand[w].min(m_solver);
            clique.add(w);
            
#ifdef CHECK_CLIQUE
            for(auto u : parallel[*v]) if(u != w and simplicial.before(*v, u)) {
                if(not parallel[u].has(w)) {
                    std::cout << "bug clique!\n";
                    
                    
                    for(auto a : relevant) {
                        std::cout << a << ":";
                        for(auto b : parallel[a]) {
                            std::cout << " " << b;
                        }
                        std::cout << std::endl;
                    }
                    
                    std::cout << "simplicial order:";
                    for (auto x{simplicial.ordering.rbegin()}; x != simplicial.ordering.rend();
                         ++x) {
                        std::cout << " " << *x;
                    }
                    std::cout << std::endl;
                    
                    
                    exit(1);
                }
            }
#endif
            
        }
    }

#ifdef DBG_CCHECK
    if (DBG_CCHECK) {
      std::cout << " = " << size_clique;
    }
#endif

    if (size_clique > capacity.max(m_solver)) {
#ifdef DBG_CCHECK
      if (DBG_CCHECK) {
        std::cout << "fail!\n";
      }
#endif
        
        bool real_clique{true};
        for (auto w{parallel[*v].begin()}; real_clique and w!=parallel[*v].end(); ++w) if(simplicial.before(*v, *w)) {
            for(auto u{w+1}; real_clique and u!=parallel[*v].end(); ++u) if(simplicial.before(*v, *u)) {
                if(not parallel[*u].has(*w)) {
                    real_clique = false;
                }
            }
        }

        if(real_clique) {
            fail_xpl.clear();
            fail_xpl.push_back(*v);
            for (auto w : parallel[*v]) {
                if (simplicial.before(*v, w))
                    fail_xpl.push_back(w);
            }
            
            
            throw Failure<T>({this, Constant::FactHint});
        } 
//        else {
//            std::cout << "do not fail because";
//            
//        }
    } else if(pruning_flag) {
        for (auto u{clique.bbegin()}; u != clique.bend(); ++u) {
            if (demand[*u].min(m_solver) + size_clique > capacity.max(m_solver)) {
                bool sequenced{false};
                bool undefined{false};
                int side{0};
                for (auto w : clique) {
                    if(end_before_start(w, *u) or end_before_start(*u, w)) {
                        sequenced = true;
                        break;
                    }
                    auto ubw{start_before_end(*u, w)};
                    auto wbu{start_before_end(w, *u)};
                    
                    side |= (ubw + 2*wbu);
                    
                    if(not ubw and not wbu) {
                        undefined = true;
                        break;
                    }
                }
                
                if(not sequenced and not undefined) {
                    if(side == 1) {
                        std::cout << prettyTask(*u) << " must be before the clique\n";
                    } else if(side == 2) {
                        std::cout << prettyTask(*u) << " must be after the clique\n";
                    } else if(side == 3) {
                        std::cout << prettyTask(*u) << " must be during the clique!\n";
                    } else {
                        std::cout << "bug\n";
                        exit(1);
                    }
                    
                    
                    for (auto w : clique) {
                        std::cout << " * " << prettyTask(w);
                        std::cout << (start_before_end(*u, w) ? " (must start before)"
                                      : (end_before_start(w, *u) ? " (must start after)" : ""));
                        std::cout << (start_before_end(w, *u) ? " (must end after)"
                                      : (end_before_start(*u, w) ? " (must end before)" : ""));
                        std::cout << std::endl;
                    }
                    std::cout << " (cap=" << capacity.max(m_solver) << ")" << std::endl;
                }
            }
        }
        
        
//      for (auto u{clique.bbegin()}; u != clique.bend(); ++u) {
//        if (demand[*u].min(m_solver) + size_clique > capacity.max(m_solver)) {
//          std::cout << prettyTask(*u) << " can't overlap with clique\n";
//          for (auto w : clique) {
//            std::cout << " * " << prettyTask(w);
//            std::cout << (start_before_end(*u, w) ? " (must start before)"
//                                                  : (end_before_start(w, *u) ? " (must start after)" : ""));
//            std::cout << (start_before_end(w, *u) ? " (must end after)"
//                                                  : (end_before_start(*u, w) ? " (must end before)" : ""));
//            std::cout << std::endl;
//
//            //                    if((start_before_end(*u, w))
//          }
//          std::cout << " (cap=" << capacity.max(m_solver) << ")" << std::endl;
//        }
//      }
    }

#ifdef DBG_CCHECK
    else if(DBG_CCHECK){
      std::cout << "\n";
    }
#endif
  }
}

//template <typename T> int CumulativeCheck<T>::getType() const {
//  return CUMULEXPL;
//}

template <typename T>
void CumulativeCheck<T>::xplain(const Literal<T> l, const hint,
                                std::vector<Literal<T>> &Cl) {
  if (l != Solver<T>::Contradiction) {
    std::cout << "bug xplain cumulative!\n";
    exit(1);
  } else {
    //        std::cout << "explain failure : clique with";
    for (auto v{fail_xpl.begin()}; v != fail_xpl.end(); ++v) {
      //            std::cout << " " << the_tasks[*v] ;
      for (auto w{v + 1}; w != fail_xpl.end(); ++w) {
        Cl.push_back(precedence[*v][*w]);
        Cl.push_back(precedence[*w][*v]);
      }
    }
    //        std::cout << std::endl;
  }
}

template <typename T>
std::ostream &CumulativeCheck<T>::display(std::ostream &os) const {
  os << "CumulativeCheck";

#ifdef DBG_CCHECK
  os << "[" << cons_id << "]";
#endif

  os << "(";
  unsigned d{0};
  for (auto &t : the_tasks) {
    std::cout << " t" << t.id() << " (" << demand[d++] << ")";
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &CumulativeCheck<T>::print_reason(std::ostream &os,
                                               const hint) const {
  //  display(os);
  os << "CumulativeCheck";
  //
  //  if (not explanations[h].empty()) {
  //
  //    auto l{explanations[h].begin()};
  //    m_solver.displayLiteral(os, *l);
  //    ++l;
  //    while (l != explanations[h].end()) {
  //      os << ", ";
  //      m_solver.displayLiteral(os, *l);
  //      ++l;
  //    }
  //  }
  //
  //  os << ")";
  return os;
}

template <typename T>
std::string CumulativeCheck<T>::prettyTask(const int i) const {
  std::stringstream ss;
  //  ss << "t" << m_tasks[i] << ": [" << est(i) << ".." << lct(i) << "] ("
  //     << minduration(i) << ")";
  ss << "t" << the_tasks[i].id() << ": ["
     << the_tasks[i].getEarliestStart(m_solver) << "..|"
    << the_tasks[i].minDuration(m_solver) << "|.."
     << the_tasks[i].getLatestEnd(m_solver) << "] (" << demand[i].min(m_solver)
     << ")";
  return ss.str();
}

// template <typename T> std::vector<int> CumulativeCheck<T>::task_map;

} // namespace tempo

#endif
