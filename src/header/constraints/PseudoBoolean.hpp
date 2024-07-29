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

#ifndef TEMPO_PseudoBooleanInterface_HPP
#define TEMPO_PseudoBooleanInterface_HPP

#include <cassert>
#include <vector>

#include "util/SparseSet.hpp"
#include "constraints/Constraint.hpp"

//#define DBG_PB

namespace tempo {

template<typename T>
class Solver;


// enforce sum(l_i) <= bound
template <typename T> class PseudoBooleanInterface : public Constraint<T> {
protected:
  Solver<T> &m_solver;

  std::vector<Literal<T>> literals;

  std::vector<T> weight;

  Reversible<T> current_bound;

  T bias;

  T ceil;

public:
  virtual T upperBound() const = 0;
  //    virtual T lowerBound() const = 0;
  virtual void setLowerBound(const T l) = 0;

  template <typename LitIter, typename WeightIter>
  PseudoBooleanInterface(Solver<T> &solver, const LitIter beg_lit,
                         const LitIter end_lit, const WeightIter beg_weight);
  virtual ~PseudoBooleanInterface();

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;
  int getType() const override;

  //    void setBound(const T ub);

  virtual std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename LitIter, typename WeightIter>
PseudoBooleanInterface<T>::PseudoBooleanInterface(Solver<T> &solver,
                                                  const LitIter beg_lit,
                                                  const LitIter end_lit,
                                                  const WeightIter beg_weight)
    : m_solver(solver), current_bound(0, &solver.getEnv()) {

  Constraint<T>::priority = Priority::High;

  bias = 0;
  ceil = 0;

  auto wp{beg_weight};
  for (auto lp{beg_lit}; lp != end_lit; ++lp) {
    auto w{*wp};
    if (w > 0) {
      literals.push_back(*lp);
      weight.push_back(w);
    } else if (w < 0) {
      literals.push_back(~(*lp));
      weight.push_back(-w);
      bias += w;
    }
    ceil += weight.back();
    ++wp;
  }

  current_bound = bias;
}

template <typename T> PseudoBooleanInterface<T>::~PseudoBooleanInterface() {}

// template <typename T> void PseudoBooleanInterface<T>::setBound(const T ub) {
//     bound = ub;
// }

template <typename T> void PseudoBooleanInterface<T>::post(const int idx) {

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

template <typename T> void PseudoBooleanInterface<T>::propagate() {

#ifdef DBG_PB
  std::cout << "beg propagate " << *this << " (lb=" << current_bound << ")"<< std::endl;
#endif

  auto i{0};
  T lb{current_bound};
  //    T lb{lowerBound()};
  for (auto p : literals) {

#ifdef DBG_PB
    if (m_solver.boolean.isUndefined(p.variable())) {
      std::cout << p << ": " << lb << " + " << weight[i] << " > "
                << upperBound() << "?\n";
    }
#endif

    if (m_solver.boolean.isUndefined(p.variable()) and
        weight[i] + lb > upperBound()) {

#ifdef DBG_PB
    std::cout << " => pruning " << ~p << std::endl;
#endif
        
      m_solver.set(~p, {this, Constant::FactHint});
    }
    ++i;
  }

#ifdef DBG_PB
  std::cout << "end propagate " << *this << std::endl;
#endif
}

template <typename T>
bool PseudoBooleanInterface<T>::notify(const Literal<T> l, const int rank) {

#ifdef DBG_PB
  std::cout << "pb notify " << l << " to " << *this << std::endl;
#endif

  if (l.isNumeric()) { // upper bound has changed, need to propagate
    return true;
  } else {
    if (m_solver.boolean.satisfied(l)) {
      // lower bound of the linear expression has changed, update the variable's
      // lb accordingly
      current_bound += weight[rank];
      setLowerBound(current_bound);

#ifdef DBG_PB
      std::cout << " -> current lb = " << current_bound << std::endl;
#endif

      return true;
    }
  }

  return false;
}

template <typename T> int PseudoBooleanInterface<T>::getType() const {
  return CARDEXPL;
}

template <typename T>
void PseudoBooleanInterface<T>::xplain(const Literal<T> l, const hint,
                                       std::vector<Literal<T>> &Cl) {
    
//    std::cout << "explain " << l << std::endl;

  auto l_lvl{(l == Solver<T>::Contradiction ? m_solver.numLiteral()
                                            : m_solver.propagationLevel(l))};
  for (auto p : literals) {
      
//      if(not m_solver.boolean.satisfied(p)) {
//          std::cout << " - "<< p << " has not contributed (not set)\n";
//      }
//      
//      if(m_solver.propagationLevel(p) >= l_lvl) {
//          std::cout << " - "<< p << " has not contributed (set later than " << l << ")\n";
//      }
      
      if (m_solver.boolean.satisfied(p) and m_solver.propagationLevel(p) < l_lvl) {
          
//          if(m_solver.propagationLevel(p) >= l_lvl) {
//              std::cout << " - "<< p << " !\n";
//          }
//          
          Cl.push_back(p);
      }
  }
}

template <typename T>
std::ostream &PseudoBooleanInterface<T>::display(std::ostream &os) const {
  os << "PB";

#ifdef DEBUG_CONSTRAINT
  os << "[" << this->id() << "]";
#endif

  os << "("; //<< weight[0] << "*" << literals[0];
  for (unsigned i{0}; i < literals.size(); ++i) {
    if (m_solver.boolean.satisfied(literals[i]))
      os << weight[i];
    else if (m_solver.boolean.falsified(literals[i]))
      os << "0";
    else
      os << weight[i] << "*" << literals[i];
    os << " + ";
  }
  os << bias << " <= ";

  //<< upperBound() << ")";
  return os;
}

template <typename T>
std::ostream &PseudoBooleanInterface<T>::print_reason(std::ostream &os,
                                                      const hint) const {
  //  display(os);
  os << "PB";
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
class PseudoBooleanConst : public PseudoBooleanInterface<T> {
public:
  template <typename lIter, typename wIter>
  PseudoBooleanConst(Solver<T> &solver, const lIter beg_lit,
                     const lIter end_lit, const wIter beg_weight, const T ub)
      : PseudoBooleanInterface<T>(solver, beg_lit, end_lit, beg_weight),
        bound(ub) {}

  T upperBound() const override { return bound; }

  //  T lowerBound() const override { return
  //  PseudoBooleanInterface<T>::current_bound; }

  void setLowerBound(const T) override{};

  virtual std::ostream &display(std::ostream &os) const override;

private:
  T bound;
  //    Reversible<T> current_bound;
};

template <typename T>
class PseudoBooleanLeqVar : public PseudoBooleanInterface<T> {
public:
  template <typename lIter, typename wIter>
  PseudoBooleanLeqVar(Solver<T> &solver, const lIter beg_lit,
                      const lIter end_lit, const wIter beg_weight,
                      const var_t ub)
      : PseudoBooleanInterface<T>(solver, beg_lit, end_lit, beg_weight),
        bound(ub) {
    //
    //    std::cout << "constr "
    //              << PseudoBooleanInterface<T>::m_solver.numeric.lower(bound)
    //              << ".."
    //              << PseudoBooleanInterface<T>::m_solver.numeric.upper(bound)
    //              << std::endl;
  }

  T upperBound() const override {
      return PseudoBooleanInterface<T>::m_solver.numeric.upper(bound);// +
//           PseudoBooleanInterface<T>::bias;
  }

  //  T lowerBound() const override {
  //    return PseudoBooleanInterface<T>::m_solver.numeric.lower(bound) +
  //    PseudoBooleanInterface<T>::bias;
  //  }

  void setLowerBound(const T l) override {
    auto p{geq<T>(bound, l)};

#ifdef DBG_PB
  std::cout << "pruning " << p << std::endl;
#endif

      PseudoBooleanInterface<T>::m_solver.set(p, {this, Constant::FactHint});
  };

  void post(const int idx) override;

  virtual std::ostream &display(std::ostream &os) const override;

private:
  var_t bound;
};

template <typename T> void PseudoBooleanLeqVar<T>::post(const int idx) {
  PseudoBooleanInterface<T>::post(idx);
  PseudoBooleanInterface<T>::m_solver.wake_me_on(ub<T>(bound), this->id());
}

template <typename T>
class PseudoBooleanGeqVar : public PseudoBooleanInterface<T> {
public:
  template <typename lIter, typename wIter>
  PseudoBooleanGeqVar(Solver<T> &solver, const lIter beg_lit,
                      const lIter end_lit, const wIter beg_weight,
                      const var_t lb)
      : PseudoBooleanInterface<T>(solver, beg_lit, end_lit, beg_weight),
        bound(lb) {
    for (auto &l : PseudoBooleanInterface<T>::literals) {
      l = ~l;
    }
    PseudoBooleanInterface<T>::ceil += PseudoBooleanInterface<T>::bias;
    PseudoBooleanInterface<T>::bias = 0;
            PseudoBooleanInterface<T>::current_bound = 0;
  }

  T upperBound() const override {
    return PseudoBooleanInterface<T>::ceil -
           PseudoBooleanInterface<T>::m_solver.numeric.lower(bound);
  }

  void setLowerBound(const T l) override {
    auto p{leq<T>(bound, PseudoBooleanInterface<T>::ceil - l)};

#ifdef DBG_PB
  std::cout << "pruning " << p << std::endl;
#endif

      PseudoBooleanInterface<T>::m_solver.set(p, {this, Constant::FactHint});
  };

  void post(const int idx) override;

  virtual std::ostream &display(std::ostream &os) const override;

private:
  var_t bound;
};

template <typename T> void PseudoBooleanGeqVar<T>::post(const int idx) {
  PseudoBooleanInterface<T>::post(idx);
  PseudoBooleanInterface<T>::m_solver.wake_me_on(lb<T>(bound), this->id());
}

template <typename T>
std::ostream &PseudoBooleanConst<T>::display(std::ostream &os) const {
  PseudoBooleanInterface<T>::display(os);
  os << bound;
  return os;
}

template <typename T>
std::ostream &PseudoBooleanLeqVar<T>::display(std::ostream &os) const {
  PseudoBooleanInterface<T>::display(os);
  os << "x" << bound << " in ["
     << PseudoBooleanInterface<T>::m_solver.numeric.lower(bound) << ".."
     << PseudoBooleanInterface<T>::m_solver.numeric.upper(bound) << "]";
  return os;
}

template <typename T>
std::ostream &PseudoBooleanGeqVar<T>::display(std::ostream &os) const {
  PseudoBooleanInterface<T>::display(os);
  os << "" << PseudoBooleanInterface<T>::ceil << " - x" << bound << " in ["
     << PseudoBooleanInterface<T>::m_solver.numeric.lower(bound) << ".."
     << PseudoBooleanInterface<T>::m_solver.numeric.upper(bound) << "]";
  return os;
}

// template <typename T> std::vector<int> PseudoBooleanInterface<T>::task_map;

} // namespace tempo

#endif
