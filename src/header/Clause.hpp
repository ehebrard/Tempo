/************************************************
 * Tempo Clause.hpp
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

#ifndef _TEMPO_CLAUSE_HPP
#define _TEMPO_CLAUSE_HPP

#include <vector>

#include "Explanation.hpp"
#include "Global.hpp"
#include "Literal.hpp"

namespace tempo {

//! Clause
/*!
 id, vector of literals and two watched
*/
template <typename T> class Clause : public std::vector<Literal<T>> {

public:
  // a constant for the empty clause
  const static Clause<T> *Empty;

  /**
   * @name constructors
   */
  //@{
  Clause(const int i = -1);
  //@}

  // a unique clause id
  int id;

  /**
   * @name watch accessors
   */
  //@{
  // the ranks of the two watched literals
  index_t watched_index[2];

  // get the two watched literals (arg in {0,1})
  Literal<T> watched(const bool) const;

  // get the rank of the literal (if it is a watched, otherwise undefined
  // behavior)
  bool watch_rank(const Literal<T>) const;

  // get the rank of the two watched literals (arg in {0,1})
  size_t watch_index(const bool) const;
  //@}

  /**
   * @name printing
   */
  //@{
  std::ostream &display(std::ostream &) const;
  //@}
};

/*!
 Implementation
*/
template <typename T>
Clause<T>::Clause(const int i) : std::vector<Literal<T>>(), id(i) {
  watched_index[0] = 0;
  watched_index[1] = 1;
}

template <typename T> Literal<T> Clause<T>::watched(const bool r) const {
  return this->operator[](watched_index[r]);
}

template <typename T> bool Clause<T>::watch_rank(const Literal<T> l) const {
  auto p{this->operator[](watched_index[1])};
  return l.sameVariable(p) and l.sign() != p.sign();
}

template <typename T> size_t Clause<T>::watch_index(const bool r) const {
  return watched_index[r];
}

template <typename T>
std::ostream &Clause<T>::display(std::ostream &os) const {
  os << id << ":";
  if (this->size() == 0)
    os << "()";
  else {
    os << "(" << this->operator[](0);
    for (size_t i{1}; i < this->size(); ++i) {
      os << ", " << this->operator[](i);
    }
      os << ")";
//      os << ") w=" << watched(0) << " & " << watched(1) << "]";
  }
  return os;
}

template <typename T>
const Clause<T> *Clause<T>::Empty = new Clause<T>(-1);

template <typename T>
std::ostream &operator<<(std::ostream &os, const Clause<T> &x) {
  return x.display(os);
}
}

#endif // _TEMPO_CLAUSE_HPP

