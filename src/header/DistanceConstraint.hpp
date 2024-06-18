#ifndef _TEMPO_DISTANCECONSTRAINT_HPP
#define _TEMPO_DISTANCECONSTRAINT_HPP

#include <iostream>

//#include "Constant.hpp"
#include "Global.hpp"

namespace tempo {

// template <typename T> class BoundSystem;

// literals x - y <= k (with x,y pointing to vars and k a constant)
template <typename T> class DistanceConstraint {

public:
  DistanceConstraint(const event f, const event t, const T d)
      : from(f), to(t), distance(d) {}

  DistanceConstraint(const var_t f, const var_t t, const T d)
      : from(static_cast<event>(f)), to(static_cast<event>(t)), distance(d) {}

  event from;
  event to;

  T distance;

  DistanceConstraint<T> operator~() const;

  static const DistanceConstraint<T> none;

  bool entails(const DistanceConstraint<T> &e) const;
  bool contradicts(const DistanceConstraint<T> &e) const;

  template <typename S> bool satisfied(S &solver) const;

  template <typename S> bool falsified(S &solver) const;

  std::ostream &display(std::ostream &os) const;
};




template <typename T>
bool operator==(const DistanceConstraint<T> &d1,
                const DistanceConstraint<T> &d2) {
  return d1.from == d2.from and d1.to == d2.to and d1.distance == d2.distance;
}

template <typename T>
const DistanceConstraint<T>
    DistanceConstraint<T>::none = DistanceConstraint<T>(-1, -1, -1);

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
  //  os << etype(to) << TASK(to) << " - " << etype(from) << TASK(from)
  //     << " <= " << distance;
  os << "x" << to << " - x" << from << " <= " << distance;
  return os;
}


template<typename T>
std::ostream &operator<<(std::ostream &os, const DistanceConstraint<T> &x) {
  return x.display(os);
}

}

#endif

