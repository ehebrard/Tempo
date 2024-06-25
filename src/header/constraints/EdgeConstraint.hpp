#ifndef TEMPO_EDGECONSTRAINT_HPP
#define TEMPO_EDGECONSTRAINT_HPP

#include <assert.h>
#include <map>
#include <vector>

#include "constraints/Constraint.hpp"
//#include "Scheduler.hpp"

//#define DBG_EDGECONS (cons_id == 0)

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

  //  int cons_id;

  //    Literal<T> r_lb; //{Constant::NoLiteral<T>};
  //    Literal<T> r_ub; //{Constant::NoLiteral<T>};
  index_t r_lb{Constant::NoIndex}; //{Constant::NoLiteral<T>};
  index_t r_ub{Constant::NoIndex}; //{Constant::NoLiteral<T>};

public:
  EdgeConstraint(Solver<T> &, const Literal<T>);
  virtual ~EdgeConstraint() = default;

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
  void xplain(const Literal<T>, const hint, std::vector<Literal<T>> &) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;
};

template <typename T>
EdgeConstraint<T>::EdgeConstraint(Solver<T> &solver, const Literal<T> p_)
    : m_solver(solver), p(p_), edge(~(m_solver.boolean.getEdge(~p)))
//,  r_lb(m_solver.numeric.strongestLiteral(bound::lower, edge.from))
//,  r_ub(m_solver.numeric.strongestLiteral(bound::upper, edge.to))
{
  // the negation of the
  Constraint<T>::priority = Priority::High;

  //
  //    std::cout << p << std::endl;
  //    std::cout << ~p << std::endl;
  //    std::cout << p.constraint() << std::endl;
  //    std::cout << (~p).constraint() << std::endl;
  //    std::cout << m_solver.boolean.getEdge(p) << std::endl;
  //    std::cout << m_solver.boolean.getEdge(~p) << std::endl;
  //    std::cout << ~(m_solver.boolean.getEdge(~p)) << std::endl;
}

template <typename T> void EdgeConstraint<T>::post(const int idx) {
  // p <=> to - from > d
  // lb(to) - ub(from) > d
  Constraint<T>::cons_id = idx;

  //    std::cout << "wake\n" << edge << std::endl;
  //    std::cout << lb<T>(edge.from) << std::endl;
  //    std::cout << ub<T>(edge.to) << std::endl;

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

    //    if (m_solver.upper(edge.to) - m_solver.lower(edge.from) <=
    //        edge.distance) {
    if (edge.satisfied(m_solver)) {

      //          r_lb = m_solver.numeric.strongestLiteral(bound::lower,
      //          edge.from); r_ub =
      //          m_solver.numeric.strongestLiteral(bound::upper, edge.to);
      r_lb = m_solver.numeric.lastLitIndex(bound::lower, edge.from);
      r_ub = m_solver.numeric.lastLitIndex(bound::upper, edge.to);

#ifdef DBG_EDGECONS
      if (DBG_EDGECONS) {
        std::cout << " ==> " << p << " (" << r_lb << "/" << r_ub << ")\n"
                  << std::endl;
      }
#endif

      m_solver.set(p, {this, NoHint});

      //            m_schedule.relax(this);
    }
  }
  //    else {
  //      m_schedule.relax(this);
  //    }

  // never propagate
  return false;
}

template <typename T> void EdgeConstraint<T>::propagate() {}

template <typename T> int EdgeConstraint<T>::getType() const {
  return EDGEEXPL;
}

template <typename T>
void EdgeConstraint<T>::xplain(const Literal<T>, const hint,
                                  std::vector<Literal<T>> &Cl) {
    
    
//    std::cout << "explain " << m_solver.pretty(p) << " by (" << edge << ") i.e, " << m_solver.getLiteral(r_lb) << " & " << m_solver.getLiteral(r_ub) << std::endl;
//    
    
//<<<<<<< HEAD
////    if(r_lb < r_ub) {
////        
////        auto p{m_solver.getLiteral(r_lb)};
//////        auto q{Literal<T>(bound::upper, edge.to, edge.distance - p.value())};
//////        auto q{geq<T>(edge.to, p.value() - edge.distance)};
////        auto q{makeNumericLiteral(bound::upper, edge.to, edge.distance - p.value())};
////        
//////        std::cout << " * or by " << p << " & " << q << std::endl << std::endl;
////        
////        Cl.push_back(p);
////        Cl.push_back(q);
////    } else {
////        
////        auto q{m_solver.getLiteral(r_ub)};
//////        auto p{Literal<T>(bound::lower, edge.from, edge.distance - q.value())};
//////        auto p{leq<T>(edge.from, edge.distance - q.value())};
////        
////        auto p{makeNumericLiteral(bound::lower, edge.from, edge.distance - q.value())};
////        
//////        std::cout << " * or by " << p << " & " << q << std::endl << std::endl;
////        
////        Cl.push_back(p);
////        Cl.push_back(q);
////    }
//=======
    if(r_lb < r_ub) {
        
        auto p{m_solver.getLiteral(r_lb)};
        auto q{makeNumericLiteral(bound::upper, edge.to, edge.distance - p.value())};
        
//        std::cout << " * or by " << p << " & " << q << std::endl << std::endl;
        
        Cl.push_back(p);
        Cl.push_back(q);
    } else {
        
        auto q{m_solver.getLiteral(r_ub)};
        auto p{makeNumericLiteral(bound::lower, edge.from, edge.distance - q.value())};
        
//        std::cout << " * or by " << p << " & " << q << std::endl << std::endl;
        
        Cl.push_back(p);
        Cl.push_back(q);
    }
//>>>>>>> 83891036076eb25d8d738038e8b2807c802515a5
//
//  Cl.push_back(m_solver.getLiteral(r_lb));
//  Cl.push_back(m_solver.getLiteral(r_ub));
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
