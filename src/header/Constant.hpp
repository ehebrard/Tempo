/************************************************
 * Tempo Constant.hpp
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

#ifndef _TEMPO_CONSTANT_HPP
#define _TEMPO_CONSTANT_HPP

#include "DistanceConstraint.hpp"
#include "Explanation.hpp"

namespace tempo {

//! Global constants
class Constant {
public:
  //  static constexpr auto DecisionHint = static_cast<hint>(0);
  static constexpr auto NoHint = static_cast<hint>(0);
  //  static constexpr auto AssumptionHint = static_cast<hint>(-1);
  //  static constexpr auto NoHint = static_cast<hint>(-1);
  static constexpr auto NoIndex = static_cast<index_t>(-1);
  static constexpr auto NoVar = static_cast<var_t>(-1);
  static constexpr index_t InfIndex = 0;
  static constexpr info_t NoSemantic = 0;
  static constexpr info_t SomeSemantic = 2;
  static constexpr var_t K = 0;
  static constexpr var_t True = 0;

  template <typename T> static constexpr T Infinity = std::numeric_limits<T>::max();

  //  template <typename T> static Explanation<T> Decision;
  template <typename T> static Explanation<T> NoReason;
  //  template <typename T> static Explanation<T> Assumption;
  template<typename T> static constexpr DistanceConstraint<T> NoEdge =
          DistanceConstraint<T>(Constant::NoVar, Constant::NoVar, Constant::Infinity<T>);

};

// template <typename T>
// Explanation<T> Constant::Decision =
//     Explanation<T>(new Explainer<T>(), Constant::DecisionHint);

template <typename T>
Explanation<T> Constant::NoReason = Explanation<T>(new Explainer<T>(),
                                                   Constant::NoHint);

// template <typename T>
// Explanation<T> Constant::Assumption =
//     Explanation<T>(new Explainer<T>(), Constant::AssumptionHint);

} // namespace tempo

#endif
