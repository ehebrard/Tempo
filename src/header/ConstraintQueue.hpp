#ifndef _TEMPO_CONSTRAINTQUEUE_HPP
#define _TEMPO_CONSTRAINTQUEUE_HPP

#include "Explanation.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"

namespace tempo {




template <typename T, int N> class ConstraintQueue {

public:
  std::vector<Constraint<T> *> &constraints;

  // which ones need to be propagated
  SparseSet<> active[N]; // alows for N levels of priority

  ConstraintQueue(std::vector<Constraint<T> *> &);
  void resize(const size_t);

  //  void initialize(std::vector<Constraint *> *c);

  // notifies constraint 'cons' of the new lit 'change' (of type "event"), @
  // position 'rank' in its scope
  void triggers(const Literal<T> l, const int rank, Constraint<T>* cons);
    void activate(const Constraint<T>* cons);

  //  // notifies constraint 'cons' of the new lit 'change' (of type "edge"), @
  //  // position 'rank' in its scope
  //  void edge_triggers(const lit var_id, const int rank, const int cons);

  // returns the active constraint of highest priority that has been
  // activated first, return NULL if there are no active constraint
  Constraint<T> *pop_front();
  // bool empty() const;

  void clear();
  bool empty() const;
  bool has(const int cons_id) const;

  std::ostream &display(std::ostream &os) const;

private:
  size_t count{0};
};

template <typename T, int N>
ConstraintQueue<T, N>::ConstraintQueue(
    std::vector<Constraint<T> *> &cons)
    : constraints(cons) {
  //    resize(constraints.size());
}

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
