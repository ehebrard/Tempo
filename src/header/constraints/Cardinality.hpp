/************************************************
 * Tempo RandomValue.hpp
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

#ifndef TEMPO_CARDINALITY_HPP
#define TEMPO_CARDINALITY_HPP

#include <cassert>
#include <vector>

#include "Solver.hpp"
#include "constraints/Constraint.hpp"

//#define DBG_LTRANS

namespace tempo {


// enforce sum(l_i) <= bound
template <typename T> class Cardinality : public Constraint<T> {
private:
  Solver<T> &m_solver;
    
  std::vector<Literal<T>> literals;

    Reversible<unsigned> current_bound;
    
    unsigned bound;

    
public:
//  template <typename Iter>
//  Cardinality(Solver<T> &solver, const Iter beg_lit,
//               const Iter end_lit, const unsigned lb);
    template <typename Iter>
    Cardinality(Solver<T> &solver, const Iter beg_var,
                 const Iter end_var, const bool sign, const unsigned lb);
  virtual ~Cardinality();

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h, std::vector<Literal<T>> &Cl) override;
  int getType() const override;
    
    void setBound(const unsigned b);

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename Iter>
Cardinality<T>::Cardinality(Solver<T> &solver, const Iter beg_var,
                            const Iter end_var, const bool sign, const unsigned b)
    : m_solver(solver), current_bound(0, &solver.getEnv()), bound(b)
{

  Constraint<T>::priority = Priority::High;

//  std::cout << "cardinality(";

  for (auto x{beg_var}; x != end_var; ++x) {
    auto l{*x == sign};

//    std::cout << " " << l;

    literals.push_back(l);
  }
//  std::cout << " ) <= " << b << std::endl;

//  current_bound = 0;
  setBound(b);
}

template <typename T> Cardinality<T>::~Cardinality() {}

template <typename T> void Cardinality<T>::setBound(const unsigned b) {
    bound = b;
}

template <typename T> void Cardinality<T>::post(const int idx) {

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
void Cardinality<T>::propagate() {}


template <typename T>
bool Cardinality<T>::notify(const Literal<T> l, const int) {
    if(m_solver.boolean.satisfied(l)) {
        ++current_bound;
        if(current_bound > bound) {
            throw Failure<T>({this, Constant::FactHint});
        } else if(current_bound == bound) {
            for(auto p : literals) {
                if(m_solver.boolean.isUndefined(p.variable()))
                    m_solver.set(~p, {this, Constant::FactHint});
            }
        }
    }
  
  return false;
}



template <typename T> int Cardinality<T>::getType() const {
  return CARDEXPL;
}

template <typename T>
void Cardinality<T>::xplain(const Literal<T> l, const hint, std::vector<Literal<T>> &Cl) {
    
    auto l_lvl{(l==Solver<T>::Contradiction ? m_solver.numLiteral() : m_solver.propagationLevel(l))};
    for(auto p : literals) {
        if(m_solver.boolean.satisfied(p) and m_solver.propagationLevel(p) < l_lvl)
            Cl.push_back(p);
    }
    
}

template <typename T>
std::ostream &Cardinality<T>::display(std::ostream &os) const {
  os << "Cardinality";

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
std::ostream &Cardinality<T>::print_reason(std::ostream &os,
                                            const hint) const {
  //  display(os);
  os << "Cardinality";
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

// template <typename T> std::vector<int> Cardinality<T>::task_map;

} // namespace tempo

#endif
