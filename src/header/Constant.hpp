
#ifndef _TEMPO_CONSTANT_HPP
#define _TEMPO_CONSTANT_HPP

#include "DistanceConstraint.hpp"
#include "Explanation.hpp"

namespace tempo {

using var_t = uint32_t;
using info_t = uint32_t;

class Constant {
public:
    
    const static index_t NoIndex;
      const static var_t NoVarx;
      const static index_t InfIndex;
    const static info_t NoSemantic;
//      const static index_t IndexOfMax;
    template <typename T> const static T Infinity;
//    template <typename T> const static T minvalue;
//    template <typename T> const static T maxvalue;
    
  static Explanation NoReason;

  template <typename T> static DistanceConstraint<T> NoEdge;
};

template <typename T>
DistanceConstraint<T> Constant::NoEdge = DistanceConstraint<T>(NOEVENT, NOEVENT, INFTY);

//template <typename T>
//const index_t Constant::no_index = static_cast<index_t>(-1);
//
//template <typename T>
//const var_t Constant::no_var = static_cast<var_t>(-1);
//
//template <typename T>
//const index_t Constant::index_of_min = 0;
//
//template <typename T>
//const index_t Constant::index_of_max = 1;

template <typename T>
const T Constant::Infinity = std::numeric_limits<T>::max();

//template <typename T>
//const T Constant::minvalue = std::numeric_limits<T>::min();

} // namespace tempo

#endif
