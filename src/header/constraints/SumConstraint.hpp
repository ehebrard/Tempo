/************************************************
 * Tempo SumConstraint.hpp
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

#ifndef TEMPO_SUMCONSTRAINT_HPP
#define TEMPO_SUMCONSTRAINT_HPP

#include <cassert>
#include <vector>

#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"

//#define DBG_SUM

namespace tempo {

template<typename T>
class Solver;

// enforce sum(w_i * x_i) <= bound
template <typename T> class SumConstraint : public Constraint<T> {
private:
  Solver<T> &m_solver;

  std::vector<var_t> scope;

  std::vector<T> weight;

  Reversible<T> current_lower_bound;

  T upper_bound;

public:
  template <typename VarIter, typename WeightIter>
  SumConstraint(Solver<T> &solver, const VarIter beg_var, const VarIter end_var,
                const WeightIter beg_weight, const T ub);
  virtual ~SumConstraint();

  bool notify(const Literal<T>, const int) override;
  void post(const int idx) override;
  void propagate() override;

  void xplain(const Literal<T> l, const hint h,
              std::vector<Literal<T>> &Cl) override;
  int getType() const override;

  std::ostream &display(std::ostream &os) const override;

  std::ostream &print_reason(std::ostream &os, const hint h) const override;

#ifdef DEBUG_CONSTRAINT
  int debug_flag{2};
#endif
};

template <typename T>
template <typename VarIter, typename WeightIter>
SumConstraint<T>::SumConstraint(Solver<T> &solver, const VarIter beg_var,
                                const VarIter end_var,
                                const WeightIter beg_weight,
                                const T ub) //, const T bias)
    : m_solver(solver), current_lower_bound(0, &solver.getEnv()),
      upper_bound(ub) {

  Constraint<T>::priority = Priority::Medium;

  //        auto wp{beg_weight};
  //  for (auto xp{beg_var}; xp != end_var; ++xp) {
  //      scope.push_back(xp->id());
  //      auto w{*wp};
  //
  //      bound -= xp->offset() * w;
  //
  //      weight.push_back(w);
  //
  //      current_lower_bound += w * (w < 0 ? m_solver.numeric.upper(scope[i]) :
  //      m_solver.numeric.lower(scope[i]));
  //
  //      ++wp;
  //  }

  //    std::cout << "new sum constraint\n";

  T lb = 0;

  auto wp{beg_weight};
  for (auto xp{beg_var}; xp != end_var; ++xp) {
    auto x{*xp};
    auto w{*wp};

    scope.push_back(x);
    weight.push_back(w);

    //        std::cout << w << "*x" << x << " in [" <<
    //        m_solver.numeric.lower(x) << ".." << m_solver.numeric.upper(x) <<
    //        "]" << std::endl;

    if (lb > -Constant::Infinity<T>) {
      if (w < 0) {
        if (m_solver.numeric.upper(x) == Constant::Infinity<T>)
          lb = -Constant::Infinity<T>;
        else
          lb += w * m_solver.numeric.upper(x);
      } else {
        if (m_solver.numeric.lower(x) == -Constant::Infinity<T>)
          lb = -Constant::Infinity<T>;
        else
          lb += w * m_solver.numeric.lower(x);
      }
    }

    ++wp;
  }

  current_lower_bound = lb;

  //    std::cout << " ==> " << static_cast<T>(current_lower_bound) <<
  //    std::endl;
}

template <typename T> SumConstraint<T>::~SumConstraint() {}

template <typename T> void SumConstraint<T>::post(const int idx) {

  Constraint<T>::cons_id = idx;
  Constraint<T>::idempotent = true;

#ifdef DEBUG_CONSTRAINT
  if (debug_flag > 0) {
    std::cout << "post " << *this << std::endl;
  }
#endif

  for (unsigned i{0}; i < scope.size(); ++i) {
    if (weight[i] < 0)
      m_solver.wake_me_on(ub<T>(scope[i]), this->id());
    else
      m_solver.wake_me_on(lb<T>(scope[i]), this->id());
  }
}

template <typename T> void SumConstraint<T>::propagate() {

#ifdef DBG_SUM
  std::cout << "propagate sum\n";
#endif

  T overall_lb{0};

  // up to one variable may have an infinite lb
  size_t inf_idx{scope.size()};

  for (unsigned i{0}; i < scope.size(); ++i) {

#ifdef DBG_SUM
    auto lb{m_solver.numeric.lower(scope[i])};
    auto ub{m_solver.numeric.upper(scope[i])};
    std::cout << scope[i] << ": " << lb
              << (lb == -Constant::Infinity<T> ? "*" : "") << ".." << ub
              << (ub == Constant::Infinity<T> ? "*" : "") << std::endl;
#endif
      
    if (weight[i] > 0) {
      if (m_solver.numeric.lower(scope[i]) == -Constant::Infinity<T>) {
        if (inf_idx < scope.size())
          return;
        inf_idx = i;
      } else {
        overall_lb += weight[i] * m_solver.numeric.lower(scope[i]);
      }
    } else {
      if (m_solver.numeric.upper(scope[i]) == Constant::Infinity<T>) {
        if (inf_idx < scope.size())
          return;
        inf_idx = i;
      } else {
        overall_lb += weight[i] * m_solver.numeric.upper(scope[i]);
      }
    }
  }

#ifdef DBG_SUM
  std::cout << weight[0] << "*x" << scope[0] << ":["
            << m_solver.numeric.lower(scope[0]) << ".."
            << m_solver.numeric.upper(scope[0]) << "]";
  for (unsigned i{1}; i < scope.size(); ++i) {
    std::cout << " + " << weight[i] << "*x" << scope[i] << ":["
              << m_solver.numeric.lower(scope[i]) << ".."
              << m_solver.numeric.upper(scope[i]) << "]";
  }
  std::cout << " <= " << upper_bound << "\nLB=" << overall_lb
            << "(inf=" << inf_idx << ")" << std::endl;
#endif

  if (inf_idx < scope.size()) {
    // only the variable with infinite lb can have its ub pruned
    auto blit{
        makeNumericLiteral((weight[inf_idx] > 0 ? bound::upper : bound::lower),
                           scope[inf_idx], overall_lb / weight[inf_idx])};

#ifdef DBG_SUM
    std::cout << " --> " << blit << std::endl;
#endif

    m_solver.set(blit, {this, Constant::FactHint});
  } else {
    for (unsigned j{0}; j < scope.size(); ++j) {
      if (weight[j] > 0) {
        T min_j{m_solver.numeric.lower(scope[j]) * weight[j]};
        T lb_except_j{overall_lb - min_j};
        T ub_j = upper_bound - lb_except_j / weight[j];

#ifdef DBG_SUM
        std::cout << " --> " << leq<T>(scope[j], ub_j) << std::endl;
#endif

        if (ub_j < m_solver.numeric.upper(scope[j])) {
          m_solver.set(leq<T>(scope[j], ub_j), {this, ub_j});
        }
      } else {
        T min_j{m_solver.numeric.upper(scope[j]) * weight[j]};

#ifdef DBG_SUM
        std::cout << " min[" << scope[j] << "] = " << min_j << std::endl;
#endif

        T lb_except_j{overall_lb - min_j};

#ifdef DBG_SUM
        std::cout << " lb\\" << scope[j] << " = " << lb_except_j << std::endl;

        std::cout << lb_except_j << " + " << weight[j] << " * lb(" << scope[j]
                  << ") <= " << upper_bound << std::endl;
#endif

        T lb_j = upper_bound - lb_except_j / weight[j];

#ifdef DBG_SUM
        std::cout << " --> " << geq<T>(scope[j], lb_j) << std::endl;
#endif

        if (lb_j > m_solver.numeric.lower(scope[j])) {
          m_solver.set(geq<T>(scope[j], lb_j), {this, lb_j});
        }
      }
    }
  }

#ifdef DBG_SUM
  std::cout << weight[0] << "*x" << scope[0] << ":["
            << m_solver.numeric.lower(scope[0]) << ".."
            << m_solver.numeric.upper(scope[0]) << "]";
  for (unsigned i{1}; i < scope.size(); ++i) {
    std::cout << " + " << weight[i] << "*x" << scope[i] << ":["
              << m_solver.numeric.lower(scope[i]) << ".."
              << m_solver.numeric.upper(scope[i]) << "]";
  }
  std::cout << std::endl;
#endif

  //    exit(1);

  //    for(unsigned j{0}; j<scope.size(); ++j) {
  //        if(weight[j] > 0) {
  //            T min_j{m_solver.numeric.lower(scope[j]) * weight[j]};
  //            T lb_except_j{overall_lb - min_j};
  //            T ub_j = upper_bound - lb_except_j/weight[j];
  //            if(ub_j < m_solver.numeric.upper(scope[j])) {
  //                m_solver.set(leq<T>(scope[j], ub_j), {this, ub_j});
  //            }
  //        } else {
  //            T min_j{m_solver.numeric.upper(scope[j]) * weight[j]};
  //            T lb_except_j{overall_lb - min_j};
  //            T lb_j = upper_bound - lb_except_j/weight[j];
  //            if(lb_j > m_solver.numeric.lower(scope[j])) {
  //                m_solver.set(geq<T>(scope[j], lb_j), {this, lb_j});
  //            }
  //        }
  //    }
}

template <typename T>
bool SumConstraint<T>::notify(const Literal<T>, const int) {

  //    std::cout << weight[0] << "*x" << scope[0] << ":[" <<
  //    m_solver.numeric.lower(scope[0])
  //    << ".." << m_solver.numeric.upper(scope[0]) << "]";
  //      for (unsigned i{1}; i<scope.size(); ++i) {
  //          std::cout << " + " << weight[i] << "*x" << scope[i] << ":[" <<
  //          m_solver.numeric.lower(scope[i])
  //          << ".." << m_solver.numeric.upper(scope[i]) << "]";
  //    }
  //    std::cout << "\nLB=" << current_lower_bound << std::endl;
  //
  ////    auto xi{scope[i]};
  //    T ai{weight[i]};
  ////    T overall_lb{static_cast<T>(current_lower_bound)};
  //
  //    auto p{m_solver.numeric.previousBound(l)};
  //
  //
  //    std::cout << "notify " << p << " -> " << l << std::endl;
  //
  //
  //    T increment = (l.value() - p.value()) * ai;
  //    if(static_cast<T>(current_lower_bound) + increment > upper_bound) {
  //        throw Failure<T>({this, Constant::FactHint});
  //    }
  //    current_lower_bound += increment;
  //
  //
  //    std::cout << "==> LB=" << current_lower_bound << std::endl;
  //
  return true;
}

template <typename T> int SumConstraint<T>::getType() const { return CARDEXPL; }

template <typename T>
void SumConstraint<T>::xplain(const Literal<T>, const hint i,
                              std::vector<Literal<T>> &Cl) {

  //    assert(scope[i] == l.variable());

  for (unsigned j{0}; j < scope.size(); ++j)
    if (static_cast<unsigned>(i) != j) {
      Cl.push_back(m_solver.numeric.getLiteral(
          (weight[i] < 0 ? bound::upper : bound::upper), scope[i]));
    }
}

template <typename T>
std::ostream &SumConstraint<T>::display(std::ostream &os) const {
  os << "SumConstraint";

#ifdef DEBUG_CONSTRAINT
  os << "[" << this->id() << "]";
#endif

  os << "(" << weight[0] << "*x" << scope[0];
  for (unsigned i{1}; i < scope.size(); ++i) {
    os << " " << weight[i] << "*x" << scope[i];
  }
  os << ")";
  return os;
}

template <typename T>
std::ostream &SumConstraint<T>::print_reason(std::ostream &os,
                                             const hint) const {
  //  display(os);
  os << "SumConstraint";
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

// template <typename T> std::vector<int> SumConstraint<T>::task_map;

} // namespace tempo

#endif
