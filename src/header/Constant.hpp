
#ifndef _TEMPO_CONSTANT_HPP
#define _TEMPO_CONSTANT_HPP

#include "DistanceConstraint.hpp"
#include "Explanation.hpp"

namespace tempo {

// template<typename T> struct Literal;

class Constant {
public:
  static constexpr auto DecisionHint = static_cast<hint>(0);
  static constexpr auto FactHint = static_cast<hint>(-1);
  static constexpr auto NoIndex = static_cast<index_t>(-1);
  static constexpr auto NoVarx = static_cast<var_t>(-1);
  static constexpr index_t InfIndex = 0;
  static constexpr info_t NoSemantic = 0;
  //      const static index_t IndexOfMax;
  template <typename T> const static T Infinity;
  //    template <typename T> const static T minvalue;
  //    template <typename T> const static T maxvalue;

//  static Explanation NoReason;
  template <typename T> static Explanation<T> Decision;
  template <typename T> static Explanation<T> GroundFact;

  template <typename T> static DistanceConstraint<T> NoEdge;

  //    template <typename T> static Literal<T> NoLiteral;
};

template <typename T>
DistanceConstraint<T> Constant::NoEdge = DistanceConstraint<T>(NOEVENT, NOEVENT, INFTY);

// template <typename T>
// Literal<T> Constant::NoLiteral = Literal<T>(Constant::NoVarx,
// std::variant<info_t,T>(static_cast<info_t>(0)));

template <typename T>
Explanation<T> Constant::Decision =
    Explanation<T>(new Explainer<T>(), Constant::DecisionHint);

template <typename T>
Explanation<T> Constant::GroundFact =
    Explanation<T>(new Explainer<T>(), Constant::FactHint);

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
const T Constant::Infinity = std::numeric_limits<T>::max(); //-Gap<T>::epsilon();

//template <typename T>
//const T Constant::minvalue = std::numeric_limits<T>::min();

} // namespace tempo

#endif
