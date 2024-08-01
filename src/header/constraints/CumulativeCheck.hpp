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
#include "Model.hpp"

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

    
    std::vector<unsigned> sorted_tasks;

public:
  template <typename ItTask, typename ItNVar, typename ItBVar>
  CumulativeCheck(Solver<T> &solver, const NumericVar<T> c, const ItTask beg_task,
               const ItTask end_task, const ItNVar beg_dem, const ItBVar beg_disj);
  virtual ~CumulativeCheck();

  bool notify(const Literal<T>, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl) override;
  int getType() const override;
    
    bool start_before_end(const int i, const int j) const;

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
                              const ItTask beg_task, const ItTask end_task,
                              const ItNVar beg_dem, ItBVar beg_disj)
    : m_solver(solver), capacity(c), relevant(std::distance(beg_task, end_task), &solver.getEnv()) {

  Constraint<T>::priority = Priority::Low;

        auto dp{beg_dem};
          for (auto jp{beg_task}; jp != end_task; ++jp) {
            the_tasks.push_back(*jp);
              demand.push_back(*dp);
              
              std::cout << *jp << " requires " << *dp << std::endl;
              ++dp;
          }

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

        std::cout << *ip << " // " << *jp << " iff (" << precedence[i][j] << " and " << precedence[j][i] << ")\n";
        
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
bool CumulativeCheck<T>::notify(const Literal<T> l, const int r) {
  
    auto x{scopex[r]};
    auto y{scopey[r]};

#ifdef DBG_CCHEK
    if (DBG_CCHEK) {
        
        //
        //    std::cout << std::endl;
        //    for (size_t i{0}; i < the_tasks.size(); ++i) {
        //      std::cout << "t" << the_tasks[i].id() << ":";
        //      for (auto j : relevant[i]) {
        //        std::cout << " t" << the_tasks[*j].id();
        //      }
        //      std::cout << std::endl;
        //    }
        
        std::cout << "\nnotify  end(" << the_tasks[x] << ") > start(" << the_tasks[y]
        << ") / " << m_solver.pretty(l) << std::endl;
        
        std::cout << start_before_end(y,x) << "|" << start_before_end(x,y) << std::endl;
        
        std::cout << "level=" << m_solver.level() << " relevant: " << relevant << std::endl;
        
        for(unsigned i{0}; i<the_tasks.size(); ++i) {
            std::cout << "t" << the_tasks[i].id()+1
            << ": ";
//            for(unsigned j{0}; j<the_tasks.size(); ++j) if(i!=j) {
//                std::cout << start_before_end(i, j);
//            } else {
//                std::cout << "*" ;
//            }
            for(unsigned j{0}; j<the_tasks.size(); ++j) 
                if(i!=j and start_before_end(i, j) and start_before_end(j, i)) {
                    std::cout << "*" ;
                    assert(parallel[i].has(j));
                } else {
                    std::cout << "." ;
                    assert(not parallel[i].has(j));
                }
            std::cout << std::endl;
        }
        
        for(unsigned i{0}; i<the_tasks.size(); ++i) {
            std::cout << "t" << the_tasks[i].id()+1
            << " ["
            << the_tasks[i].getEarliestStart(m_solver)
            << "-" << the_tasks[i].getLatestStart(m_solver)
            << ".."
            << the_tasks[i].getEarliestEnd(m_solver)
            << "-" << the_tasks[i].getLatestEnd(m_solver)
            << "] p=" << the_tasks[i].minDuration(m_solver)
            << " d=" <<
            demand[i].min(m_solver)
            << std::endl;
        }
    }
#endif
    
//    std::cout << m_solver.pretty(precedence[x][y]) << " <> "
//    << m_solver.pretty(precedence[y][x]) << std::endl;
    
//    bool new_relevant{false};
//    if(m_solver.boolean.satisfied(precedence[y][x])) {
    if(start_before_end(x,y)) {
        if(not relevant.has(x)) {
            relevant.add(x);
//            new_relevant = true;
        }
        if(not relevant.has(y)) {
            relevant.add(y);
//            new_relevant = true;
        }
        
        assert(not parallel[x].has(y) and not parallel[y].has(x));
        parallel[x].add(y);
        parallel[y].add(x);
        
        return true;
//        std::cout << "relevant: " << relevant << std::endl;
    }
    
#ifdef DBG_CCHEK
    if (DBG_CCHEK) {
        std::cout << "no new relevantism\n";
    }
#endif
    
    
  return false;
}


template <typename T> void CumulativeCheck<T>::propagate() {
    
#ifdef DBG_CCHEK
    if (DBG_CCHEK) {
        std::cout << "propagate " <<  *this << std::endl;
        std::cout << "relevant: " << relevant << std::endl;
    }
#endif
    
    assert(not relevant.empty());
    
    sorted_tasks.clear();
    //    for(unsigned i{0}; i<the_tasks.size(); ++i)
    //        sorted_tasks.push_back(i);
    for(auto i : relevant) {
        sorted_tasks.push_back(i);
    }
    
    
    //    std::sort(sorted_tasks.begin(), sorted_tasks.end(),
    //              [&](const int x, const int y) -> bool {
    //        return the_tasks[x].getEarliestEnd(m_solver) <
    //        the_tasks[y].getEarliestEnd(m_solver);
    //    });
    
    std::sort(sorted_tasks.begin(), sorted_tasks.end(),
              [&](const int x, const int y) -> bool {
        return not start_before_end(x,y);
    });
    
#ifdef DBG_CCHEK
    if (DBG_CCHEK) {
        for(auto i : sorted_tasks) {
            std::cout << "t" << the_tasks[i].id()+1
            << " ["
            << the_tasks[i].getEarliestStart(m_solver)
            << "-" << the_tasks[i].getLatestStart(m_solver)
            << ".."
            << the_tasks[i].getEarliestEnd(m_solver)
            << "-" << the_tasks[i].getLatestEnd(m_solver)
            << "] p=" << the_tasks[i].minDuration(m_solver)
            << " d=" <<
            //            m_solver.numeric.lower(demand[i])
            demand[i].min(m_solver)
//            << ": ";
//            for(auto j : sorted_tasks) if(i!=j) {
//                if(start_before_end(i, j)) {
//                    if(start_before_end(j, i))
//                        std::cout << "*";
//                    else
//                        std::cout << ">";
//                } else {
//                    if(start_before_end(j, i))
//                        std::cout << "<";
//                    else
//                        std::cout << ".";
//                }
//            }
//            std::cout 
            << std::endl;
        }
        
        for(unsigned i{0}; i<the_tasks.size(); ++i) {
            std::cout << "t" << the_tasks[i].id()+1 << ": ";
            for(unsigned j{0}; j<the_tasks.size(); ++j)
                if(i!=j and start_before_end(i, j) and start_before_end(j, i)) {
                    std::cout << "*" ;
                    assert(parallel[i].has(j));
                } else {
                    std::cout << "." ;
                    assert(not parallel[i].has(j));
                }
            std::cout << std::endl;
        }
        
        for(unsigned i{0}; i<the_tasks.size(); ++i) {
            std::cout << "t" << the_tasks[i].id()+1 << ": ";
            for(unsigned j{0}; j<the_tasks.size(); ++j)
                if(i!=j) {
                    std::cout << start_before_end(i, j) ;
                } else {
                    std::cout << "." ;
                }
            std::cout << std::endl;
        }
    }
    
    std::vector<index_t> neighbors;
#endif
    
    
    
    while(not sorted_tasks.empty()) {
        auto i{sorted_tasks.back()};
        sorted_tasks.pop_back();
        
        T size_clique{demand[i].min(m_solver)};
        
#ifdef DBG_CCHEK
        if (DBG_CCHEK) {
            std::cout << "clique: " << the_tasks[i].id()+1 << " (" << size_clique << ")";
        }
#endif
        
        
        for(auto j : sorted_tasks) {
            if(start_before_end(i,j) and start_before_end(j,i)) {
                size_clique += demand[j].min(m_solver);
                
#ifdef DBG_CCHEK
                if (DBG_CCHEK) {
                    
                    std::cout << " + " << the_tasks[j].id()+1 << " (" << size_clique << ")";
                    
                    for(auto k : neighbors) if(j != k)
                        if(not start_before_end(j,k) or not start_before_end(k,j)) {
                            std::cout << std::endl << the_tasks[i].id()+1 << " is not simplicial because "
                            << the_tasks[j].id()+1 << " and " <<  the_tasks[k].id()+1 << " are not adjacent\n";
                            exit(1);
                        }
                    neighbors.push_back(j);
                }
#endif
                
                if(size_clique > capacity.max(m_solver)) {
                    
#ifdef DBG_CCHEK
                if (DBG_CCHEK) {
                    std::cout << "fail!\n";
                }
#endif
                    
                    throw Failure<T>({this,Constant::FactHint});
                }
            }
        }
        
        neighbors.clear();
#ifdef DBG_CCHEK
                if (DBG_CCHEK) {
                    std::cout << "\n";
                }
#endif
        
    }
    
    
    
    
    
//    auto start_clique{sorted_tasks.begin()};
//    auto end_clique{start_clique};
//    T size_clique{demand[*start_clique].min(m_solver)};
//    while(++end_clique != sorted_tasks.end()) {
//        
//#ifdef DBG_CCHEK
//        if (DBG_CCHEK) {
//            std::cout << "\nadd t" << the_tasks[*end_clique].id()+1 << "\n";
//        }
//#endif
//        
//        size_clique += demand[*end_clique].min(m_solver);
//        
//#ifdef DBG_CCHEK
//        if (DBG_CCHEK) {
//            std::cout << "size=" << size_clique << " need to reduce clique " << "?\n";
//            
//            if(start_before_end(*end_clique, *start_clique))
//                std::cout << "t" << the_tasks[*end_clique].id()+1 << " must start before "
//                << "t" << the_tasks[*start_clique].id()+1 << " ends" << std::endl;
//            else
//                std::cout << "t" << the_tasks[*end_clique].id()+1 << " may start after "
//                << "t" << the_tasks[*start_clique].id()+1 << " ends" << std::endl;
//        }
//#endif
//        
//        while(not (start_before_end(*end_clique, *start_clique) and start_before_end(*start_clique, *end_clique))) {
//            size_clique -= demand[*start_clique].min(m_solver);
//#ifdef DBG_CCHEK
//            if (DBG_CCHEK) {
//                std::cout << "remove t" << the_tasks[*start_clique].id()+1 << ", size=" << size_clique << "\n";
//            }
//#endif
//            if(++start_clique == end_clique) {
//                break;
//            }
//            
//#ifdef DBG_CCHEK
//            if (DBG_CCHEK) {
//                if(start_before_end(*end_clique, *start_clique))
//                    std::cout << "t" << the_tasks[*end_clique].id()+1 << " must start before "
//                    << "t" << the_tasks[*start_clique].id()+1 << " ends" << std::endl;
//                else
//                    std::cout << "t" << the_tasks[*end_clique].id()+1 << " may start after "
//                    << "t" << the_tasks[*start_clique].id()+1 << " ends" << std::endl;
//                
//                if(start_before_end(*start_clique, *end_clique))
//                    std::cout << "t" << the_tasks[*start_clique].id()+1 << " must start before "
//                    << "t" << the_tasks[*end_clique].id()+1 << " ends" << std::endl;
//                else
//                    std::cout << "t" << the_tasks[*start_clique].id()+1 << " may start after "
//                    << "t" << the_tasks[*end_clique].id()+1 << " ends" << std::endl;
//            }
//#endif
//        }
//        
//        
//#ifdef DBG_CCHEK
//        if (DBG_CCHEK) {
//            auto i{start_clique};
//            
//            std::cout << "\ncurrent clique (size=" << size_clique << "):\n";
//            do {
//                std::cout << "t" << the_tasks[*i].id()+1 << " [" << the_tasks[*i].getEarliestStart(m_solver) << ".."
//                << the_tasks[*i].getLatestEnd(m_solver) << "] p=" << the_tasks[*i].minDuration(m_solver)
//                << " d=" <<
//                //            m_solver.numeric.lower(demand[*i])
//                demand[*i].min(m_solver)
//                << std::endl;
//            } while(i++ != end_clique);
//        }
//#endif
//        
//        //        if(size_clique > m_solver.numeric.upper(capacity))
//        if(size_clique > capacity.max(m_solver)) {
//            
//#ifdef DBG_CCHEK
//            if (DBG_CCHEK) {
//                std::cout << "FAIL!\n";
//            }
//#endif
//            
//            throw Failure<T>({this,Constant::FactHint});
//        }
//    }
}

template <typename T> int CumulativeCheck<T>::getType() const {
  return CUMULEXPL;
}

template <typename T>
void CumulativeCheck<T>::xplain(const Literal<T>, const hint,
                             std::vector<Literal<T>> &) {
}

template <typename T>
std::ostream &CumulativeCheck<T>::display(std::ostream &os) const {
  os << "CumulativeCheck";

#ifdef DEBUG_CONSTRAINT
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
     << the_tasks[i].getEarliestStart(m_solver) << ".."
     << the_tasks[i].getLatestEnd(m_solver) << "] ("
     << the_tasks[i].minDuration(m_solver) << ")";
  return ss.str();
}

// template <typename T> std::vector<int> CumulativeCheck<T>::task_map;

} // namespace tempo

#endif
