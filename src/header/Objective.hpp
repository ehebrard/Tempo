/************************************************
 * Tempo Objective.hpp
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

#ifndef _TEMPO_OBJECTIVE_HPP
#define _TEMPO_OBJECTIVE_HPP

namespace tempo {

template <typename T> class Solver;

//! \class Objective
/*!
 \brief Wraper for objective
 Objective is either minimization or maximization of a given numeric variable
 In order to model fancy objective, add a variable; use constraints to define
 the objective as this variable; and minimize or maximize it
*/
template <typename T> class Objective {
public:
  Objective(NumericVar<T> &x) : X(x) {}

  T gap() { return p_b - d_b; }
  T dualBound() const { return d_b; }
  T primalBound() const { return p_b; }

  void setDual(const T v) { d_b = v; }

protected:
  NumericVar<T> &X;
  T d_b{0};
  T p_b{Constant::Infinity<T>};
};

template <typename T> class MinimizationObjective : public Objective<T> {
public:
  MinimizationObjective(NumericVar<T> &x) : Objective<T>(x) {}

  T value(Solver<T> &solver) { return Objective<T>::X.min(solver); }

  void setPrimal(const T v, Solver<T> &solver) {
    Objective<T>::p_b = v;
    if (Objective<T>::gap()) {
      apply(Objective<T>::p_b, solver);
    }
  }

private:
  void apply(const T target, Solver<T> &solver) {
    //    solver.set(Objective<T>::X < target);
    solver.set(Objective<T>::X.before(target - Gap<T>::epsilon()));
  }
};

template <typename T> class MaximizationObjective : public Objective<T> {
public:
  MaximizationObjective(NumericVar<T> &x) : Objective<T>(x) {}

  T value(Solver<T> &solver) { return Objective<T>::X.max(solver); }

  void setPrimal(const T v, Solver<T> &solver) {
    Objective<T>::p_b = v;
    if (Objective<T>::gap()) {
      apply(Objective<T>::p_b, solver);
    }
  }

private:
  void apply(const T target, Solver<T> &solver) {
    solver.set(Objective<T>::X.after(target + Gap<T>::epsilon()));
  }
};

} // namespace tempo

#endif
