/************************************************
 * Tempo OptionalEdgeConstraint.hpp
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

#ifndef TEMPO_OPTIONALEDGECONSTRAINT_HPP
#define TEMPO_OPTIONALEDGECONSTRAINT_HPP

#include <assert.h>
#include <map>
#include <vector>

#include "constraints/Constraint.hpp"

namespace tempo {

template <typename T> class Solver;
template <typename T> class DistanceConstraint;


// constraint that watches events e_i and e_j and sets the corresponding edge to
// false
template <typename T> class OptionalEdgeConstraint : public Constraint<T> {
private:
  Solver<T> &m_solver;

  // the literal to propagate for each edge
  Literal<T> p[2];
  // the edges themselves
  DistanceConstraint<T> edge[2];
    // constraint edge[false] implies literal p[false] and
    // constraint edge[true] implies literal p[true]
    
    var_t exist;

  index_t r_lb[2] = {Constant::NoIndex,Constant::NoIndex};
    index_t r_ub[2] = {Constant::NoIndex,Constant::NoIndex};

public:
    OptionalEdgeConstraint(Solver<T> &, const var_t, const var_t);
  virtual ~OptionalEdgeConstraint() = default;

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;
};

template <typename T>
OptionalEdgeConstraint<T>::OptionalEdgeConstraint(Solver<T> &solver, const var_t x, const var_t e)
: m_solver(solver)
{
  Constraint<T>::priority = Priority::High;
    
    p[false] = m_solver.boolean.getLiteral(true, x);
    p[true] = m_solver.boolean.getLiteral(false, x);
    edge[false] = ~(m_solver.boolean.getEdge(p[true]));
    edge[true] = ~(m_solver.boolean.getEdge(p[false]));
    exist = e;
}

template <typename T> void OptionalEdgeConstraint<T>::post(const int idx) {
  // p <=> to - from > d
  // lb(to) - ub(from) > d
  Constraint<T>::cons_id = idx;
  m_solver.wake_me_on(lb<T>(edge[false].from), this->id());
  m_solver.wake_me_on(ub<T>(edge[false].to), this->id());
    m_solver.wake_me_on(lb<T>(edge[true].from), this->id());
    m_solver.wake_me_on(ub<T>(edge[true].to), this->id());
    if(exist != Constant::NoVar)
        m_solver.wake_me_on(makeBooleanLiteral<T>(true,exist,0), this->id());
}

template <typename T>
bool OptionalEdgeConstraint<T>::notify(const Literal<T>, const int rank) {

#ifdef DBG_EDGECONS
  if (DBG_EDGECONS) {
    std::cout << *this << "\nnotified of " << l << ": "
              << "x" << edge[false].to << " <= " << m_solver.numeric.upper(edge[false].to)
              << " & x" << edge[false].from
              << " >= " << m_solver.numeric.lower(edge[false].from)
              << " <= " << edge[false].distance
      << " AND x" << edge[true].to << " <= " << m_solver.numeric.upper(edge[true].to)
      << " & x" << edge[true].from
      << " >= " << m_solver.numeric.lower(edge[true].from)
      << " <= " << edge[true].distance << "\n";
  }
#endif
  
    if(rank < 4) {
        bool i = (rank >= 2);
        
        if (m_solver.boolean.isUndefined(p[i].variable())) {
            
            
            auto relevant{exist == Constant::NoVar or m_solver.boolean.isTrue(exist)};
            if (edge[i].satisfied(m_solver)) {
                
                if(relevant) {
                    r_lb[i] = m_solver.numeric.lastLitIndex(bound::lower, edge[i].from);
                    r_ub[i] = m_solver.numeric.lastLitIndex(bound::upper, edge[i].to);
                    
#ifdef DBG_EDGECONS
                    if (DBG_EDGECONS) {
                        std::cout << " ==> " << p[i] << " (" << r_lb[i] << "/" << r_ub[i] << ")\n"
                        << std::endl;
                    }
#endif
                    
                    m_solver.set(p[i], {this, i});
                } else if(edge[1-i].satisfied(m_solver)) {
                    for(auto i{0}; i<2; ++i) {
                        r_lb[i] = m_solver.numeric.lastLitIndex(bound::lower, edge[i].from);
                        r_ub[i] = m_solver.numeric.lastLitIndex(bound::upper, edge[i].to);
                    }
                    m_solver.set(makeBooleanLiteral<T>(false, exist, 0), {this, i});
                }
            }
        }
    } else {
        for(auto i{0}; i<2; ++i) {
            if (edge[i].satisfied(m_solver)) {
                r_lb[i] = m_solver.numeric.lastLitIndex(bound::lower, edge[i].from);
                r_ub[i] = m_solver.numeric.lastLitIndex(bound::upper, edge[i].to);
                m_solver.set(p[i], {this, i});
            }
        }
    }

  // never propagate
  return false;
}

template <typename T> void OptionalEdgeConstraint<T>::propagate() {}

template <typename T> int OptionalEdgeConstraint<T>::getType() const {
  return EDGEEXPL;
}

template <typename T>
void OptionalEdgeConstraint<T>::xplain(const Literal<T> l, const hint h,
                                  std::vector<Literal<T>> &Cl) {
    
    if(l.hasSemantic()) {
        
        if(r_lb[h] < r_ub[h]) {
            
            auto a{m_solver.getLiteral(r_lb[h])};
            auto b{makeNumericLiteral(bound::upper, edge[h].to, edge[h].distance - a.value())};
            
            Cl.push_back(a);
            Cl.push_back(b);
            if(exist != Constant::NoVar)
                Cl.push_back(makeBooleanLiteral<T>(true, exist, 0));
            
        } else {
            
            auto a{m_solver.getLiteral(r_ub[h])};
            auto b{makeNumericLiteral(bound::lower, edge[h].from, edge[h].distance - a.value())};
            
            Cl.push_back(a);
            Cl.push_back(b);
            if(exist != Constant::NoVar)
                Cl.push_back(makeBooleanLiteral<T>(true, exist, 0));
        }
    } else {
        // explain removal
        for(auto i{0}; i<2; ++i) {
            Literal<T> a;
            Literal<T> b;
            if(r_lb[i] < r_ub[i]) {
                a = m_solver.getLiteral(r_lb[i]);
                b = makeNumericLiteral(bound::upper, edge[i].to, edge[i].distance - a.value());
            } else {
                a = m_solver.getLiteral(r_ub[h]);
                b = makeNumericLiteral(bound::lower, edge[h].from, edge[h].distance - a.value());
            }
            Cl.push_back(a);
            Cl.push_back(b);
        }
        
    }

}

template <typename T>
std::ostream &OptionalEdgeConstraint<T>::display(std::ostream &os) const {
  os << "watcher [" << edge[false] << "] => " << p[false] << " and [" << edge[true] << "] => " << p[true] ;
  return os;
}

template <typename T>
std::ostream &OptionalEdgeConstraint<T>::print_reason(std::ostream &os,
                                                 const hint) const {
  os << "[" << edge << "]";
  return os;
}

} // namespace tempo

#endif
