#ifndef TEMPO_EDGECONSTRAINT_HPP
#define TEMPO_EDGECONSTRAINT_HPP

#include <assert.h>
#include <map>
#include <vector>

#include "constraints/Constraint.hpp"
//#include "Scheduler.hpp"

//#define DBG_EDGECONS (cons_id == 0)

namespace tempo {

template <typename T> class Scheduler;
template <typename T> class DistanceConstraint;

// constraint that watches events e_i and e_j and sets the corresponding edge to
// false
template <typename T> class EdgeConstraint : public Constraint {
private:
  Scheduler<T> &m_schedule;

  // the literal to propagate for the edge
  const lit p;
  // the edge itself
  const DistanceConstraint<T> edge;

  //  int cons_id;

  genlit r_lb{NoLit};
  genlit r_ub{NoLit};

public:
  EdgeConstraint(Scheduler<T> &, const lit);
  virtual ~EdgeConstraint() = default;

  bool notify_bound(const lit, const int) override;
  void post(const int idx) override;
  void propagate() override;

  // void xplain(const lit l, const hint h, std::vector<lit> &Cl) const
  // override;
  void xplain(const lit, const hint, std::vector<lit> &) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;
};

template <typename T>
EdgeConstraint<T>::EdgeConstraint(Scheduler<T> &scheduler, const lit p_)
    : m_schedule(scheduler), p(p_), edge(~(m_schedule.getEdge(NOT(p)))) {
  // the negation of the
  priority = Priority::High;
}

// template <typename T>
// void EdgeConstraint<T>::post(const int idx) {
//     // (from, to, d) <=> to - from <= d
//     //
//     // ub(to) - lb(from) <= d
//     cons_id = idx;
//     m_schedule.wake_me_on_event(LOWERBOUND(edge.from), cons_id);
//     m_schedule.wake_me_on_event(UPPERBOUND(edge.to), cons_id);
// }
//
// template <typename T>
// bool EdgeConstraint<T>::notify_bound(const lit, const int) {
//
//     if(m_schedule.upper(edge.to) - m_schedule.lower(edge.from) <=
//     edge.distance)
//         m_schedule.set(p, {this,
//         static_cast<hint>(m_schedule.numBoundLiteral())});
//
//     // never propagate
//     return false;
// }

template <typename T> void EdgeConstraint<T>::post(const int idx) {
  // p <=> to - from > d
  // lb(to) - ub(from) > d
  Constraint::cons_id = idx;
  m_schedule.wake_me_on_event(LOWERBOUND(edge.from), cons_id);
  m_schedule.wake_me_on_event(UPPERBOUND(edge.to), cons_id);
}

template <typename T>
bool EdgeConstraint<T>::notify_bound(const lit, const int) {

#ifdef DBG_EDGECONS
  if (DBG_EDGECONS) {
    std::cout << *this << "\nnotified of " << prettyEventLit(l) << ": "
              << prettyEvent(edge.to) << " <= " << m_schedule.upper(edge.to)
              << " & " << prettyEvent(edge.from)
              << " >= " << m_schedule.lower(edge.from)
              << " <= " << edge.distance << "\n";
  }
#endif

  if (m_schedule.isUndefined(VAR(p))) {

    if (m_schedule.upper(edge.to) - m_schedule.lower(edge.from) <=
        edge.distance) {

      r_lb = m_schedule.getBoundIndex(LOWERBOUND(edge.from));
      r_ub = m_schedule.getBoundIndex(UPPERBOUND(edge.to));

#ifdef DBG_EDGECONS
      if (DBG_EDGECONS) {
        std::cout << " ==> " << m_schedule.prettyLiteral(EDGE(p)) << " ("
                  << r_lb << "/" << r_ub << ")\n";
        std::cout << " #lit = " << m_schedule.numBoundLiteral() << std::endl;
      }
#endif

      m_schedule.set(p, {this, NoHint});

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
void EdgeConstraint<T>::xplain(const lit, const hint, std::vector<lit> &Cl) {
  ////    if(l == NoLit) {
  ////        "?";
  ////    } else {
  ////        "?";
  ////    }
  //    std::cout << "TODO!";

  //    assert(l != NoLit);
  //    assert(LTYPE(l) == EDGE_LIT);
  //    assert(FROM_GEN(l) == p);

  //    std::cout << cons_id << ": " << r_lb << " / " << r_ub << std::endl;
  //    std::cout << " #lit = " << m_schedule.numBoundLiteral() << std::endl;

  Cl.push_back(BOUND(r_lb));
  Cl.push_back(BOUND(r_ub));
}

template <typename T>
std::ostream &EdgeConstraint<T>::display(std::ostream &os) const {
  os << "watcher [" << edge << "] => [" << m_schedule.getEdge(p) << "]";
  return os;
}

template <typename T>
std::ostream &EdgeConstraint<T>::print_reason(std::ostream &os,
                                              const hint) const {
  os << "not[" << edge << "]";
  return os;
}

} // namespace tempo

#endif
