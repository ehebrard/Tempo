
#ifndef _TEMPO_CONSTANT_HPP
#define _TEMPO_CONSTANT_HPP

#include "DistanceConstraint.hpp"
#include "Explanation.hpp"

namespace tempo {

class Constant {
public:
  static Explanation NoReason;

  template <typename T> static DistanceConstraint<T> NoEdge;
};

template <typename T>
DistanceConstraint<T> Constant::NoEdge = DistanceConstraint<T>(NOEVENT, NOEVENT,
                                                               INFTY);

} // namespace tempo

#endif
