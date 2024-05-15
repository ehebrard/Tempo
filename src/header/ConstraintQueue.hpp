#ifndef _TEMPO_CONSTRAINTQUEUE_HPP
#define _TEMPO_CONSTRAINTQUEUE_HPP

#include "Explanation.hpp"
#include "constraints/Constraint.hpp"
#include "util/SparseSet.hpp"

namespace tempo {


template <int N> class ConstraintQueue {

public:
  std::vector<Constraint *>& constraints;

  // which ones need to be propagated
  SparseSet<> active[N]; // alows for N levels of priority

  ConstraintQueue(std::vector<Constraint*>&);
    void resize(const size_t);

//  void initialize(std::vector<Constraint *> *c);

  // notifies constraint 'cons' of the new lit 'change' (of type "event"), @ position 'rank' in its scope
    void bound_triggers(const lit var_id, const int rank, const int cons,
                        const Explainer *responsible);

    // notifies constraint 'cons' of the new lit 'change' (of type "edge"), @ position 'rank' in its scope
    void edge_triggers(const lit var_id, const int rank, const int cons,
                       const Explainer *responsible);

    // returns the active constraint of highest priority that has been
    // activated first, return NULL if there are no active constraint
    Constraint *pop_front();
    // bool empty() const;

    void clear();
    bool empty() const;
    bool has(const int cons_id) const;

    std::ostream &display(std::ostream &os) const;

  private:
    size_t count{0};
};

template <int N>
ConstraintQueue<N>::ConstraintQueue(std::vector<Constraint*>& cons) : constraints(cons) {
//    resize(constraints.size());
}

template <int N>
void ConstraintQueue<N>::resize(const size_t n) {
    for (auto i{0}; i < N; ++i)
        active[i].reserve(n);
}


//template <int N>
//void ConstraintQueue<N>::initialize(std::vector<Constraint*> *c) {
//  constraints = c;
//  for (auto i{0}; i < N; ++i)
//    active[i].reserve(constraints->size());
//}

// notifies that variable 'id' has changed, activate the corresponding
// triggers
template <int N>
void ConstraintQueue<N>::bound_triggers(const lit var_id, const int rank,
                                        const int cons_id,
                                        const Explainer *responsible) {

  //  assert(cons_id >= 0);
  //  assert(cons_id < static_cast<int>(constraints->size()));
  //    assert(active[cons->priority].size() == constraints->size());

  auto cons = constraints[cons_id];

  //    if(cons->idempotent)
  //        std::cout << *responsible << std::endl;

  //  if (responsible == static_cast<Explainer *>(cons))
  //    std::cout << "it happens!\n";
  //
  //  if (responsible->id() == cons->id()) {
  //    std::cout << "here!\n";
  //  }

  if (not cons->idempotent or responsible != static_cast<Explainer *>(cons)) {
    if (cons->notify_bound(var_id, rank) and
        not active[to_underlying(cons->priority)].has(cons_id)) {
      active[to_underlying(cons->priority)].add(cons_id);
      ++count;
    }
  }
}

// notifies that variable 'id' has changed, activate the corresponding
// triggers
template <int N>
void ConstraintQueue<N>::edge_triggers(const lit var_id, const int rank,
                                       const int cons_id,
                                       const Explainer *responsible) {

  assert(cons_id >= 0);
  assert(cons_id < static_cast<int>(constraints.size()));

  auto cons = constraints[cons_id];

  //  if (responsible == static_cast<Explainer *>(cons))
  //    std::cout << "it happens!\n";
  //
  //  if (responsible->id() == cons->id()) {
  //    std::cout << "here!\n";
  //  }

  if (not cons->idempotent or responsible != static_cast<Explainer *>(cons)) {

    if (cons->notify_edge(var_id, rank) and
        not active[to_underlying(cons->priority)].has(cons_id)) {
      active[to_underlying(cons->priority)].add(cons_id);
      ++count;
    }
  }
}

template <int N>
Constraint *ConstraintQueue<N>::pop_front() {
  auto i{N};
  Constraint *cons{NULL};
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

template <int N>
void ConstraintQueue<N>::clear() {
  for (auto i{N}; i-- > 0;) {
    active[i].clear();
  }
    count = 0;
}

template <int N> bool ConstraintQueue<N>::empty() const {
//  for (auto i{N}; i-- > 0;)
//    if (not active[i].empty())
//      return false;
//  return true;
    return count==0;
}

template <int N>
bool ConstraintQueue<N>::has(const int cons_id) const {
//  for (auto i{N}; i-- > 0;)
//    if (active[i].has(cons_id))
//      return true;
//  return false;
    return active[constraints[cons_id]->priority].has(cons_id);
}

template <int N>
std::ostream &ConstraintQueue<N>::display(std::ostream &os) const {
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

template <int N>
std::ostream &operator<<(std::ostream &os, const ConstraintQueue<N> &x) {
  return x.display(os);
}

template <typename T, int N> class NewConstraintQueue {

public:
  std::vector<NewConstraint<T> *> &constraints;

  // which ones need to be propagated
  SparseSet<> active[N]; // alows for N levels of priority

  NewConstraintQueue(std::vector<NewConstraint<T> *> &);
  void resize(const size_t);

  //  void initialize(std::vector<Constraint *> *c);

  // notifies constraint 'cons' of the new lit 'change' (of type "event"), @
  // position 'rank' in its scope
  void bound_triggers(const lit var_id, const int rank, const int cons);

  // notifies constraint 'cons' of the new lit 'change' (of type "edge"), @
  // position 'rank' in its scope
  void edge_triggers(const lit var_id, const int rank, const int cons);

  // returns the active constraint of highest priority that has been
  // activated first, return NULL if there are no active constraint
  NewConstraint<T> *pop_front();
  // bool empty() const;

  void clear();
  bool empty() const;
  bool has(const int cons_id) const;

  std::ostream &display(std::ostream &os) const;

private:
  size_t count{0};
};

template <typename T, int N>
NewConstraintQueue<T, N>::NewConstraintQueue(
    std::vector<NewConstraint<T> *> &cons)
    : constraints(cons) {
  //    resize(constraints.size());
}

template <typename T, int N>
void NewConstraintQueue<T, N>::resize(const size_t n) {
  for (auto i{0}; i < N; ++i)
    active[i].reserve(n);
}

// template <typename T, int N>
// void NewConstraintQueue<T,N>::initialize(std::vector<Constraint*> *c) {
//   constraints = c;
//   for (auto i{0}; i < N; ++i)
//     active[i].reserve(constraints->size());
// }

// notifies that variable 'id' has changed, activate the corresponding
// triggers
template <typename T, int N>
void NewConstraintQueue<T, N>::bound_triggers(const lit var_id, const int rank,
                                              const int cons_id) {

  //  assert(cons_id >= 0);
  //  assert(cons_id < static_cast<int>(constraints->size()));
  //    assert(active[cons->priority].size() == constraints->size());

  auto cons = constraints[cons_id];

  //    if(cons->idempotent)
  //        std::cout << *responsible << std::endl;

  //  if (responsible == static_cast<Explainer *>(cons))
  //    std::cout << "it happens!\n";
  //
  //  if (responsible->id() == cons->id()) {
  //    std::cout << "here!\n";
  //  }

  //  if (not cons->idempotent or responsible != static_cast<Explainer *>(cons))
  //  {
  if (cons->notify_bound(var_id, rank) and
      not active[to_underlying(cons->priority)].has(cons_id)) {
    active[to_underlying(cons->priority)].add(cons_id);
    ++count;
  }
  //  }
}

// notifies that variable 'id' has changed, activate the corresponding
// triggers
template <typename T, int N>
void NewConstraintQueue<T, N>::edge_triggers(const lit var_id, const int rank,
                                             const int cons_id) {

  assert(cons_id >= 0);
  assert(cons_id < static_cast<int>(constraints.size()));

  auto cons = constraints[cons_id];

  //  if (responsible == static_cast<Explainer *>(cons))
  //    std::cout << "it happens!\n";
  //
  //  if (responsible->id() == cons->id()) {
  //    std::cout << "here!\n";
  //  }

  //  if (not cons->idempotent or responsible != static_cast<Explainer *>(cons))
  //  {

  if (cons->notify_edge(var_id, rank) and
      not active[to_underlying(cons->priority)].has(cons_id)) {
    active[to_underlying(cons->priority)].add(cons_id);
    ++count;
  }
  //  }
}

template <typename T, int N>
NewConstraint<T> *NewConstraintQueue<T, N>::pop_front() {
  auto i{N};
  NewConstraint<T> *cons{NULL};
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

template <typename T, int N> void NewConstraintQueue<T, N>::clear() {
  for (auto i{N}; i-- > 0;) {
    active[i].clear();
  }
  count = 0;
}

template <typename T, int N> bool NewConstraintQueue<T, N>::empty() const {
  //  for (auto i{N}; i-- > 0;)
  //    if (not active[i].empty())
  //      return false;
  //  return true;
  return count == 0;
}

template <typename T, int N>
bool NewConstraintQueue<T, N>::has(const int cons_id) const {
  //  for (auto i{N}; i-- > 0;)
  //    if (active[i].has(cons_id))
  //      return true;
  //  return false;
  return active[constraints[cons_id]->priority].has(cons_id);
}

template <typename T, int N>
std::ostream &NewConstraintQueue<T, N>::display(std::ostream &os) const {
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
std::ostream &operator<<(std::ostream &os, const NewConstraintQueue<T, N> &x) {
  return x.display(os);
}

} // namespace tempo

#endif // _TEMPO_CONSTRAINTQUEUE_HPP
