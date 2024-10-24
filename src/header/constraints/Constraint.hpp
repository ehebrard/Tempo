/************************************************
 * Tempo Constraint.hpp
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

#ifndef _TEMPO_CONSTRAINT_HPP
#define _TEMPO_CONSTRAINT_HPP

#include <ostream>

#include "Explanation.hpp"
#include "Global.hpp"

namespace tempo {


template <typename T> class Constraint : public Explainer<T> {

public:

    // priority on the constraint queue
  Priority priority = Priority::High;
    
    // whether the constraint should be called on literals it is responsible for
  bool idempotent{false};

  // give the constraint the explainer id 'idx' (this is where it should subscribe to literals)
  virtual void post(const int idx) = 0;
  // propagate the constraint
  virtual void propagate() = 0;
  // notify a change (with the literal and it's variable rank in the scope)
  virtual bool notify(const Literal<T>, const int) { return false; }

  virtual std::ostream &display(std::ostream &os) const = 0;

};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Constraint<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif // _TEMPO_CONSTRAINT_HPP
