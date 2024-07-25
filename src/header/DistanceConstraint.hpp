/************************************************
 * Tempo DistanceConstraint.hpp
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

#ifndef _TEMPO_DISTANCECONSTRAINT_HPP
#define _TEMPO_DISTANCECONSTRAINT_HPP

#include <iostream>

#include "Global.hpp"

namespace tempo {

//! Difference logic constraint
/*!
 literals  x - y \leq k (with x,y pointing to vars and k a constant)
  - the corresponding edge in the temporal graph is y --> x with label k
 */
template <typename T> class DistanceConstraint {

public:
  DistanceConstraint() {}

  DistanceConstraint(const var_t f, const var_t t, const T d)
      : from(f), to(t), distance(d) {}

  var_t from;
  var_t to;

  T distance;

  DistanceConstraint<T> operator~() const;

//  static const DistanceConstraint<T> none;

  bool entails(const DistanceConstraint<T> &e) const;
  bool contradicts(const DistanceConstraint<T> &e) const;
    
    bool isNull() const { return from == static_cast<var_t>(-1); }

  template <typename S> bool satisfied(S &solver) const;

  template <typename S> bool falsified(S &solver) const;

  std::ostream &display(std::ostream &os) const;
};




template <typename T>
bool operator==(const DistanceConstraint<T> &d1,
                const DistanceConstraint<T> &d2) {
  return d1.from == d2.from and d1.to == d2.to and d1.distance == d2.distance;
}

//template <typename T>
//const DistanceConstraint<T>
//    DistanceConstraint<T>::none = DistanceConstraint<T>(-1, -1, -1);

template <typename T>
DistanceConstraint<T> DistanceConstraint<T>::operator~() const {
  return {to, from, -distance - Gap<T>::epsilon()};
}

template <typename T>
template <typename S>
bool DistanceConstraint<T>::satisfied(S &solver) const {
  return solver.numeric.upper(to) - solver.numeric.lower(from) <= distance;
}

template <typename T>
template <typename S>
bool DistanceConstraint<T>::falsified(S &solver) const {
  return solver.numeric.lower(to) - solver.numeric.upper(from) > distance;
}

template <typename T>
bool DistanceConstraint<T>::entails(const DistanceConstraint<T> &e) const {
  return e.from == from and e.to == to and distance <= e.distance;
}

template <typename T>
bool DistanceConstraint<T>::contradicts(const DistanceConstraint<T> &e) const {
  return e.from == to and e.to == from and e.distance + distance < 0;
}

template <typename T>
std::ostream &DistanceConstraint<T>::display(std::ostream &os) const {
  os << "x" << to << " - x" << from << " <= " << distance;
  return os;
}


template<typename T>
std::ostream &operator<<(std::ostream &os, const DistanceConstraint<T> &x) {
  return x.display(os);
}

}

#endif

