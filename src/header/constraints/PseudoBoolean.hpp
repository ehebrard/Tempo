/************************************************
 * Tempo PseudoBoolean.hpp
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

#ifndef TEMPO_PSEUDOBOOLEAN_HPP
#define TEMPO_PSEUDOBOOLEAN_HPP

#include <cassert>
#include <vector>

#include "Solver.hpp"
#include "util/SparseSet.hpp"
#include "constraints/Constraint.hpp"

//#define DBG_LTRANS

namespace tempo {


// enforce sum(l_i) <= bound
template <typename T> class PseudoBoolean : public Constraint<T> {
private:
  Solver<T> &m_solver;
    
  std::vector<Literal<T>> literals;
    
    std::vector<T> weight;

    Reversible<T> current_bound;
    
    T bound;
    
    T bias;
    
//    SparseSet<> scope;

    
public:
    template <typename LitIter, typename WeightIter>
    PseudoBoolean(Solver<T> &solver, const LitIter beg_lit,
                  const LitIter end_lit, const WeightIter beg_weight, const T ub);//, const T bias=0);
  virtual ~PseudoBoolean();

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl) override;
  int getType() const override;
    
    void setBound(const T ub);

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename LitIter, typename WeightIter>
PseudoBoolean<T>::PseudoBoolean(Solver<T> &solver, const LitIter beg_lit,
             const LitIter end_lit, const WeightIter beg_weight, const T ub)//, const T bias)
: m_solver(solver), current_bound(0, &solver.getEnv()), bound(ub), bias(0)
{

  Constraint<T>::priority = Priority::High;
        
        auto wp{beg_weight};
  for (auto lp{beg_lit}; lp != end_lit; ++lp) {
      auto w{*wp};
      if(w>0) {
          literals.push_back(*lp);
          weight.push_back(w);
      } else if(w<0) {
          literals.push_back(~(*lp));
          weight.push_back(-(w));
          bias += w;
      }
      ++wp;
  }
        
    current_bound = bias;
}

template <typename T> PseudoBoolean<T>::~PseudoBoolean() {}

template <typename T> void PseudoBoolean<T>::setBound(const T ub) {
    bound = ub;
}

template <typename T> void PseudoBoolean<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
    Constraint<T>::idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

    for(auto l : literals) {
        m_solver.wake_me_on(l, this->id());
    }
    
}

template <typename T>
void PseudoBoolean<T>::propagate() {
    auto i{0};
    T lb{current_bound};
    for(auto p : literals) {
        if(m_solver.boolean.isUndefined(p.variable()) and weight[i] + lb > bound)
            m_solver.set(~p, {this, Constant::FactHint});
    }
}


template <typename T>
bool PseudoBoolean<T>::notify(const Literal<T> l, const int rank) {
    if(m_solver.boolean.satisfied(l)) {
        current_bound += weight[rank];
        if(current_bound > bound) {
            throw Failure<T>({this, Constant::FactHint});
        } else {
            return true;
        }
    }
  
  return false;
}



template <typename T> int PseudoBoolean<T>::getType() const {
  return CARDEXPL;
}

template <typename T>
void PseudoBoolean<T>::xplain(const Literal<T> l, const hint, std::vector<Literal<T>> &Cl) {
    
    auto l_lvl{(l==Solver<T>::Contradiction ? m_solver.numLiteral() : m_solver.propagationLevel(l))};
    for(auto p : literals) {
        if(m_solver.boolean.satisfied(p) and m_solver.propagationLevel(p) < l_lvl)
            Cl.push_back(p);
    }
    
}

template <typename T>
std::ostream &PseudoBoolean<T>::display(std::ostream &os) const {
  os << "PseudoBoolean";

#ifdef DEBUG_CONSTRAINT
  os << "[" << this->id() << "]";
#endif

  os << "(" << literals[0];
    for (unsigned i{1}; i<literals.size(); ++i) {
      std::cout << " " << literals[i];
  }
  std::cout << ")";
  return os;
}

template <typename T>
std::ostream &PseudoBoolean<T>::print_reason(std::ostream &os,
                                            const hint) const {
  //  display(os);
  os << "PseudoBoolean";
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

// template <typename T> std::vector<int> PseudoBoolean<T>::task_map;

} // namespace tempo

#endif
