/************************************************
 * Tempo EdgeConstraint.hpp
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

#ifndef TEMPO_EDGECONSTRAINT_HPP
#define TEMPO_EDGECONSTRAINT_HPP

#include <assert.h>
#include <map>
#include <vector>

#include "constraints/Constraint.hpp"

namespace tempo {

template <typename T> class Solver;
template <typename T> class DistanceConstraint;


// constraint that watches events e_i and e_j and sets the corresponding edge to
// false
template <typename T> class EdgeConstraint : public Constraint<T> {
private:
  Solver<T> &m_solver;

  // the literal to propagate for the edge
  const Literal<T> p;
  // the edge itself
  const DistanceConstraint<T> edge;

  index_t r_lb{Constant::NoIndex};
  index_t r_ub{Constant::NoIndex};

public:
  EdgeConstraint(Solver<T> &, const Literal<T>);
  virtual ~EdgeConstraint() = default;

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
//  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;
};

template <typename T>
EdgeConstraint<T>::EdgeConstraint(Solver<T> &solver, const Literal<T> p_)
    : m_solver(solver), p(p_), edge(~(m_solver.boolean.getEdge(~p)))
{
  Constraint<T>::priority = Priority::High;
}

template <typename T> void EdgeConstraint<T>::post(const int idx) {
  // p <=> to - from > d
  // lb(to) - ub(from) > d
  Constraint<T>::cons_id = idx;
  m_solver.wake_me_on(lb<T>(edge.from), Constraint<T>::cons_id);
  m_solver.wake_me_on(ub<T>(edge.to), Constraint<T>::cons_id);
}

template <typename T>
bool EdgeConstraint<T>::notify(const Literal<T>, const int) {

#ifdef DBG_EDGECONS
  if (DBG_EDGECONS) {
    std::cout << *this << "\nnotified of " << l << ": "
              << "x" << edge.to << " <= " << m_solver.numeric.upper(edge.to)
              << " & x" << edge.from
              << " >= " << m_solver.numeric.lower(edge.from)
              << " <= " << edge.distance << "\n";
  }
#endif

  if (m_solver.boolean.isUndefined(p.variable())) {

    if (edge.satisfied(m_solver)) {
        
      r_lb = m_solver.numeric.lastLitIndex(bound::lower, edge.from);
      r_ub = m_solver.numeric.lastLitIndex(bound::upper, edge.to);

#ifdef DBG_EDGECONS
      if (DBG_EDGECONS) {
        std::cout << " ==> " << p << " (" << r_lb << "/" << r_ub << ")\n"
                  << std::endl;
      }
#endif

      m_solver.set(p, {this, Constant::FactHint});
    }
  }
//      else {
//        m_solver.relax(this);
//          //optimal    4299017    8993174    36868     17    459206    18   243.928
//          //optimal    3884559    8125410    35409     17    680596    18   229.468
//      }

  // never propagate
  return false;
}

template <typename T> void EdgeConstraint<T>::propagate() {}

//template <typename T> int EdgeConstraint<T>::getType() const {
//  return EDGEEXPL;
//}

template <typename T>
void EdgeConstraint<T>::xplain(const Literal<T>, const hint,
                                  std::vector<Literal<T>> &Cl) {
    
    if(r_lb < r_ub) {
        
        auto p{m_solver.getLiteral(r_lb)};
        auto q{makeNumericLiteral(bound::upper, edge.to, edge.distance - p.value())};
         
        Cl.push_back(p);
        Cl.push_back(q);
    } else {
        
        auto q{m_solver.getLiteral(r_ub)};
        auto p{makeNumericLiteral(bound::lower, edge.from, edge.distance - q.value())};
         
        Cl.push_back(p);
        Cl.push_back(q);
    }

}

template <typename T>
std::ostream &EdgeConstraint<T>::display(std::ostream &os) const {
  os << "watcher [" << edge << "] => [" << p << "]";
  return os;
}

template <typename T>
std::ostream &EdgeConstraint<T>::print_reason(std::ostream &os,
                                                 const hint) const {
  os << "[" << edge << "]";
  return os;
}

} // namespace tempo

#endif
