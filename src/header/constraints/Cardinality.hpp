/************************************************
 * Tempo CardinalityInterface.hpp
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

#include "constraints/Constraint.hpp"
#include "ReversibleObject.hpp"

//#define DBG_CARD

namespace tempo {

template<typename T>
class Solver;


// enforce sum(l_i) <= bound
template <typename T> class CardinalityInterface : public Constraint<T> {
protected:
  Solver<T> &m_solver;
    
  std::vector<Literal<T>> literals;

    Reversible<T> current_bound;

    
public:
    
    virtual T upperBound() = 0;
    virtual T lowerBound() = 0;
    virtual void setLowerBound(const T l) = 0;
    
//  template <typename Iter>
//  CardinalityInterface(Solver<T> &solver, const Iter beg_lit,
//               const Iter end_lit, const unsigned lb);
  template <typename Iter>
  CardinalityInterface(Solver<T> &solver, const Iter beg_var, const Iter end_var,
              const bool sign);

  template <typename Iter>
  CardinalityInterface(Solver<T> &solver, const Iter beg_lit, const Iter end_lit);
  virtual ~CardinalityInterface();

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;
//  int getType() const override;

//  void setBound(const unsigned b);

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename Iter>
CardinalityInterface<T>::CardinalityInterface(Solver<T> &solver, const Iter beg_var,
                            const Iter end_var, const bool sign)
    : m_solver(solver), current_bound(0, &solver.getEnv())
{

  Constraint<T>::priority = Priority::High;
  for (auto x{beg_var}; x != end_var; ++x) {
    auto l{*x == sign};
    literals.push_back(l);
  }
//  setBound(b);
}

template <typename T>
template <typename Iter>
CardinalityInterface<T>::CardinalityInterface(Solver<T> &solver, const Iter beg_lit,
                            const Iter end_lit)
    : m_solver(solver), current_bound(0, &solver.getEnv()) {

  Constraint<T>::priority = Priority::High;
  for (auto l{beg_lit}; l != end_lit; ++l) {
    literals.push_back(*l);
  }
//  setBound(b);
}

template <typename T> CardinalityInterface<T>::~CardinalityInterface() {}

//template <typename T> void CardinalityInterface<T>::setBound(const unsigned b) {
//    bound = b;
//}

template <typename T> void CardinalityInterface<T>::post(const int idx) {

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
void CardinalityInterface<T>::propagate() {}


template <typename T>
bool CardinalityInterface<T>::notify(const Literal<T> l, const int) {
    
    auto ub{upperBound()};

    if(l.isNumeric()) {
        if(current_bound > ub) {
            throw Failure<T>({this, Constant::FactHint});
        } else if(current_bound == ub) {
            for(auto p : literals) {
                if(m_solver.boolean.isUndefined(p.variable()))
                    m_solver.set(~p, {this, Constant::FactHint});
            }
        }
    } else {
        
        if(m_solver.boolean.satisfied(l)) {
            ++current_bound;

            T lb{static_cast<T>(current_bound)};
            if(lb > ub) {
                throw Failure<T>({this, Constant::FactHint});
            } else {
              if (lb > lowerBound()) {
                setLowerBound(lb);
              }
              if (lb == ub) {
                for (auto p : literals) {
                  if (m_solver.boolean.isUndefined(p.variable()))
                    m_solver.set(~p, {this, Constant::FactHint});
                }
              }
            }
        }
    }

  return false;
}



//template <typename T> int CardinalityInterface<T>::getType() const {
//  return CARDEXPL;
//}

template <typename T>
void CardinalityInterface<T>::xplain(const Literal<T> l, const hint, std::vector<Literal<T>> &Cl) {
    
    auto l_lvl{(l==Solver<T>::Contradiction ? m_solver.numLiteral() : m_solver.propagationLevel(l))};
    for(auto p : literals) {
        if(m_solver.boolean.satisfied(p) and m_solver.propagationLevel(p) < l_lvl)
            Cl.push_back(p);
    }
    
}

template <typename T>
std::ostream &CardinalityInterface<T>::display(std::ostream &os) const {
  os << "Cardinality";

#ifdef DEBUG_CONSTRAINT
  os << "[" << this->id() << "]";
#endif

    os << "(" ;
    if(m_solver.boolean.satisfied(literals[0]))
        std::cout << "1";
    else if(m_solver.boolean.falsified(literals[0]))
        std::cout << "0";
    else
        std::cout << literals[0];
    for (unsigned i{1}; i<literals.size(); ++i) {
        if(m_solver.boolean.satisfied(literals[i]))
            std::cout << " 1";
        else if(m_solver.boolean.falsified(literals[i]))
            std::cout << " 0";
        else
            std::cout << " " << literals[i];
  }
  std::cout << ")";
  return os;
}

template <typename T>
std::ostream &CardinalityInterface<T>::print_reason(std::ostream &os,
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



template <typename T> class CardinalityConst : public CardinalityInterface<T> {
public:
    
    template <typename Iter>
    CardinalityConst(Solver<T> &solver, const Iter beg_var, const Iter end_var,
                const bool sign, const T ub) : CardinalityInterface<T>(solver, beg_var, end_var, sign), bound(ub) {}

    template <typename Iter>
    CardinalityConst(Solver<T> &solver, const Iter beg_lit, const Iter end_lit,
                const T ub) : CardinalityInterface<T>(solver, beg_lit, end_lit), bound(ub) {}
    
    T upperBound() override {return bound;}
    
    T lowerBound() override {return CardinalityInterface<T>::current_bound;}
    
    void setLowerBound(const T) override {};
    
private:
    T bound;
};


template <typename T> class CardinalityLeqVar : public CardinalityInterface<T> {
public:
    
    template <typename Iter>
    CardinalityLeqVar(Solver<T> &solver, const Iter beg_var, const Iter end_var,
                const bool sign, const var_t ub) : CardinalityInterface<T>(solver, beg_var, end_var, sign), bound(ub) {}

    template <typename Iter>
    CardinalityLeqVar(Solver<T> &solver, const Iter beg_lit, const Iter end_lit,
                const var_t ub) : CardinalityInterface<T>(solver, beg_lit, end_lit), bound(ub) {}
    
    T upperBound() override {return CardinalityInterface<T>::m_solver.numeric.upper(bound);}
    
    T lowerBound() override {return CardinalityInterface<T>::m_solver.numeric.lower(bound);}
//    T lowerBound() override {return CardinalityInterface<T>::current_bound;}
    
    void setLowerBound(const T l) override {
        CardinalityInterface<T>::m_solver.set(geq<T>(bound, l));
//        
//        std::cout << " ==> set [" << CardinalityInterface<T>::m_solver.numeric.lower(bound)
//        << ".." << CardinalityInterface<T>::m_solver.numeric.upper(bound) << "]\n";
    };
    
    void post(const int idx) override;
    
//    void propagate() override;
    
private:
    var_t bound;
};


template <typename T> void CardinalityLeqVar<T>::post(const int idx) {
    CardinalityInterface<T>::post(idx);
    CardinalityInterface<T>::m_solver.wake_me_on(ub<T>(bound), this->id());
}



template <typename T> class CardinalityGeqVar : public CardinalityInterface<T> {
public:
    
    template <typename Iter>
    CardinalityGeqVar(Solver<T> &solver, const Iter beg_var, const Iter end_var,
                const bool sign, const var_t lb) : CardinalityInterface<T>(solver, beg_var, end_var, sign), bound(lb) {}

    template <typename Iter>
    CardinalityGeqVar(Solver<T> &solver, const Iter beg_lit, const Iter end_lit,
                const var_t lb) : CardinalityInterface<T>(solver, beg_lit, end_lit), bound(lb) {
        for(auto& l : CardinalityInterface<T>::literals) {
            l = ~l;
        }
    }
    
    T upperBound() override {return static_cast<T>(CardinalityInterface<T>::literals.size()) -  CardinalityInterface<T>::m_solver.numeric.lower(bound);}
//        T upperBound() override {return CardinalityInterface<T>::m_solver.numeric.upper(bound);}
    
//    T lowerBound() override {return CardinalityInterface<T>::current_bound;}
    T lowerBound() override {return static_cast<T>(CardinalityInterface<T>::literals.size()) - CardinalityInterface<T>::m_solver.numeric.upper(bound);}
    
    void setLowerBound(const T l) override {
        
//        std::cout << " set [" << CardinalityInterface<T>::m_solver.numeric.lower(bound)
//        << ".." << CardinalityInterface<T>::m_solver.numeric.upper(bound) << "] <= " 
//        << (l - static_cast<T>(CardinalityInterface<T>::literals.size())) << "\n";
//        
        CardinalityInterface<T>::m_solver.set(leq<T>(bound, static_cast<T>(CardinalityInterface<T>::literals.size()) - l));
//
//            std::cout << " ==> [" << CardinalityInterface<T>::m_solver.numeric.lower(bound)
//            << ".." << CardinalityInterface<T>::m_solver.numeric.upper(bound) << "] (" << lowerBound() << ")\n";
    };
    
    void post(const int idx) override;
    
private:
    var_t bound;
};


template <typename T> void CardinalityGeqVar<T>::post(const int idx) {
    CardinalityInterface<T>::post(idx);
    CardinalityInterface<T>::m_solver.wake_me_on(lb<T>(bound), this->id());
}

// template <typename T> std::vector<int> CardinalityInterface<T>::task_map;

} // namespace tempo

#endif
