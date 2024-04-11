#ifndef TEMPO_CARDINALITY_HPP
#define TEMPO_CARDINALITY_HPP

#include <cassert>
#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "ReversibleObject.hpp"
#include "Scheduler.hpp"
#include "constraints/Constraint.hpp"

//#define DBG_LTRANS

namespace tempo {

template <typename T> class CardinalityGeq : public Constraint {
private:
  Scheduler<T> &m_schedule;

  std::vector<lit> literals;

    Reversible<size_t> current_ub;
    
    size_t lower_limit;
    
    std::vector<lit> expl;
    
public:
  template <typename Iter>
  CardinalityGeq(Scheduler<T> &scheduler, const Iter beg_lit,
               const Iter end_lit, const size_t lb);
  virtual ~CardinalityGeq();

  bool notify_edge(const int lit, const int rank) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const lit l, const hint h, std::vector<lit> &Cl) override;
  int getType() const override;
    
    void setBound(const int l);

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename Iter>
CardinalityGeq<T>::CardinalityGeq(Scheduler<T> &scheduler, const Iter beg_lit,
                            const Iter end_lit, const size_t lb)
    : m_schedule(scheduler), lower_limit(lb)
{

  priority = LOW;

  for (auto l{beg_lit}; l != end_lit; ++l) {
      literals.push_back(*l);
  }
        
        current_ub = literals.size();
        setBound(lb);
}

template <typename T> CardinalityGeq<T>::~CardinalityGeq() {}

template <typename T> void CardinalityGeq<T>::setBound(const int l) {
    lower_limit = l;
}

template <typename T> void CardinalityGeq<T>::post(const int idx) {

  cons_id = idx;
  idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

    for(auto l : literals) {
        m_schedule.wake_me_on_edge(l, cons_id);
    }
    
}

template <typename T>
void CardinalityGeq<T>::propagate() {}


template <typename T>
bool CardinalityGeq<T>::notify_edge(const lit l, const int) {
    if(m_schedule.falsified(l)) {
        --current_ub;
        if(current_ub < lower_limit) {
            expl.clear();
            for(auto p : literals) {
                if(m_schedule.isFalse(VAR(p)))
                    expl.push_back(EDGE(p));
            }
            assert(expl.size() > static_cast<size_t>(lower_limit));
            throw Failure({this, NoHint}
                          
                          );
        } else if(current_ub == lower_limit) {
            expl.clear();
            for(auto p : literals) {
                if(m_schedule.isUndefined(VAR(p)))
                    m_schedule.set(p);
                else if(m_schedule.isFalse(VAR(p))) {
                    expl.push_back(EDGE(p));
                }
            }
            assert(expl.size() == static_cast<size_t>(lower_limit));
        }
    }
  
  return false;
}



template <typename T> int CardinalityGeq<T>::getType() const {
  return CARDEXPL;
}

template <typename T>
void CardinalityGeq<T>::xplain(const lit, const hint, std::vector<lit> &Cl) {
    
    for(auto l : expl) {
        Cl.push_back(l);
    }
    
}

template <typename T>
std::ostream &CardinalityGeq<T>::display(std::ostream &os) const {
  os << "CardinalityGeq";

#ifdef DEBUG_CONSTRAINT
  os << "[" << cons_id << "]";
#endif

  os << "(";
  for (auto l : literals) {
      std::cout << " " << m_schedule.prettyLiteral(l);
  }
  std::cout << " )";
  return os;
}

template <typename T>
std::ostream &CardinalityGeq<T>::print_reason(std::ostream &os,
                                            const hint) const {
  //  display(os);
  os << "CardinalityGeq";
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

// template <typename T> std::vector<int> CardinalityGeq<T>::task_map;

} // namespace tempo

#endif
