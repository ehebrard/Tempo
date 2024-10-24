/************************************************
 * Tempo Solver.hpp
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

#ifndef _TEMPO_CONSTRAINTQUEUE_HPP
#define _TEMPO_CONSTRAINTQUEUE_HPP

#include "Explanation.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"
#include "util/enum.hpp"

namespace tempo {

//! Prioirity queue for constraints
/*!
Constraints are triggered ['notify()'] on literal they subscribed to
 - if they have been successfuly triggered once ('notify()' returned true), they
are enqueued
 - they are dequeued on demand in an order compatible with their priority
declaration
 - there are N priority levels (N is a parameter)
*/
template <typename T, int N> class ConstraintQueue {

public:
  /**
   * @name constructors
   */
  //@{
  ConstraintQueue(std::vector<Constraint<T> *> &);
  void resize(const size_t);
  //@}

  /**
   * @name accessors
   */
  //@{
  // notifies constraint 'cons' of the new lit 'change' (of type "event"), @
  // position 'rank' in its scope
  void triggers(const Literal<T> l, const int rank, Constraint<T>* cons);
  // enqueue a constraint (if it is not yet in the queue)
  void activate(const Constraint<T> *cons);

  // returns the active constraint of highest priority that has been
  // activated first, return NULL if there are no active constraint
  Constraint<T> *pop_front();

  // clear the queue
  void clear();

  // whether there the queue is empty
  bool empty() const;

  // whether the queue contains the constraint with id 'cons_id'
  bool has(const int cons_id) const;
  //@}

  /**
   * @name printing
   */
  //@{
  std::ostream &display(std::ostream &os) const;
  //@}

private:
  // pointer to the constraint list
  std::vector<Constraint<T> *> &constraints;

  // which ones need to be propagated
  SparseSet<> active[N]; // alows for N levels of priority

  // to speed-up the call to 'empty()' when there are many priority levels
  size_t count{0};
};

/*!
 Implementation
*/
template <typename T, int N>
ConstraintQueue<T, N>::ConstraintQueue(std::vector<Constraint<T> *> &cons)
    : constraints(cons) {}

template <typename T, int N>
void ConstraintQueue<T, N>::resize(const size_t n) {
  for (auto i{0}; i < N; ++i)
    active[i].reserve(n);
}

// template <typename T, int N>
// void ConstraintQueue<T,N>::initialize(std::vector<Constraint*> *c) {
//   constraints = c;
//   for (auto i{0}; i < N; ++i)
//     active[i].reserve(constraints->size());
// }

// notifies that variable 'id' has changed, activate the corresponding
// triggers
template <typename T, int N>
void ConstraintQueue<T, N>::triggers(const Literal<T> l, const int rank,
                                         Constraint<T>* cons) {

//  auto cons = constraints[cons_id];

//  if (cons->notify(l, rank) and
//      not active[to_underlying(cons->priority)].has(cons_id)) {
//    active[to_underlying(cons->priority)].add(cons_id);
//    ++count;
//  }
    if (cons->notify(l, rank)) {
        activate(cons);
    }
    
}

template <typename T, int N>
void ConstraintQueue<T, N>::activate(const Constraint<T>* cons) {
    auto cons_id{cons->id()};
  if (not active[to_underlying(cons->priority)].has(cons_id)) {
    active[to_underlying(cons->priority)].add(cons_id);
    ++count;
  }
}

template <typename T, int N>
Constraint<T> *ConstraintQueue<T, N>::pop_front() {
  auto i{N};
  Constraint<T> *cons{NULL};
  while (i-- > 0) {
    if (not active[i].empty()) {
      cons = constraints[active[i].front()];
      active[i].pop_front();
      --count;
      break;
    }
  }
  return cons;
}

template <typename T, int N> void ConstraintQueue<T, N>::clear() {
  for (auto i{N}; i-- > 0;) {
    active[i].clear();
  }
  count = 0;
}

template <typename T, int N> bool ConstraintQueue<T, N>::empty() const {
  //  for (auto i{N}; i-- > 0;)
  //    if (not active[i].empty())
  //      return false;
  //  return true;
  return count == 0;
}

template <typename T, int N>
bool ConstraintQueue<T, N>::has(const int cons_id) const {
  //  for (auto i{N}; i-- > 0;)
  //    if (active[i].has(cons_id))
  //      return true;
  //  return false;
  return active[constraints[cons_id]->priority].has(cons_id);
}

template <typename T, int N>
std::ostream &ConstraintQueue<T, N>::display(std::ostream &os) const {
  size_t check{0};
  for (auto p{N}; p-- > 0;) {
    os << "p" << p << ":";
    for (auto a : active[p])
      os << " " << *(constraints[a]);
    os << std::endl;
    check += active[p].size();
  }
  assert(count == check);
  return os;
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, const ConstraintQueue<T, N> &x) {
  return x.display(os);
}

} // namespace tempo

#endif // _TEMPO_CONSTRAINTQUEUE_HPP
