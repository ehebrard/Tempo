/************************************************
 * Tempo Global.hpp
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

#ifndef _TEMPO_GLOBAL_HPP
#define _TEMPO_GLOBAL_HPP

#include <string>
#include <limits>
#include <cstdint>
#include <sys/resource.h>

//! Global definitions
namespace tempo {

// primitive type for variables (variables are essentially indices)
using var_t = std::uint32_t;
// primitive type for the datafield of literals (type, sign, variable)
using info_t = std::uint32_t;
// primitive type for indexing anything except variables
using index_t = std::uint32_t;
// for Boolean variables that can be undefined [used only for the output of solver.satisfiable()?]
using boolean_state = int;
#define TrueState 1
#define FalseState 0
#define UnknownState -1
// primitive type given to an explanation algorithm together with the literal to explain [used to encode the relevant info to be able to explain]
using hint = int;
//#define NoHint -1
//
//#define DBG_BOUND num_choicepoints >= 0
//#define DBG_CBOUND solver.num_choicepoints >= 0
//#define DBG_BBOUND num_choicepoints >= 0 //(sched.num_fails >= 236)
//#define DBG_TRACE 1                      // 183 //1+2+4+32+128
#define SEARCH 1
#define DOMAINS 2
#define BRANCH 4
#define CLAUSES 8
#define PROPAGATION 16
#define LEARNING 32
#define QUEUE 64
#define UNITPROPAGATION 128

//#define DBG_MINIMIZATION
//#define DBG_EDGEFINDING (m_solver.num_cons_propagations >= 0)
//#define DBG_EXPLEF true      //(m_schedule.num_fails > 236)
//#define DBG_THETA (m_schedule.num_fails >= 481)
//#define DBG_BELLMAN true //(sched.num_choicepoints >= 1045)
//#define DBG_BELLMAN_EXPL (sched.num_choicepoints >= 4172)
//#define DEBUG_HEURISTICS
//#define DBG_UP
//#define DBG_CL 10000000
//#define
//#define DBG_CLPLUS 1030
//#define DBG_TRANSITIVITY true //(m_schedule.num_choicepoints >= 4064)
//#define DBG_EXPL_TRANS true
//#define DBG_SOL
//#define DBG_FAIL true
//#define DBG_CCHECK m_solver.num_choicepoints >= 8620
//#define DBG_LEXBFS true
//#define DBG_FTRANS true //m_solver.num_choicepoints >= 900
//#define DBG_BELLMAN_FT true
//#define DBG_EXPL_FTRANS true
//#define DBG_SEF false //(this->id() == 321)
//#define DBG_EXTRACT true

// priority values for constraint propagation
enum class Priority {
    Low = 0,
    Medium,
    High
};

/**
 * Converts enum to underlying type
 * @tparam E enum type
 * @param e enum to convert
 * @return value of underlying type
 * @note replace with std impl in c++23: https://en.cppreference.com/w/cpp/utility/to_underlying
 */
template<typename E>
constexpr auto to_underlying(E e) noexcept {
    return static_cast<std::underlying_type_t<E>>(e);
}

//#define FTRANSEXPL 10
//#define CUMULEXPL 9
//#define CARDEXPL 8
//#define TRANSITIVITYEXPL 7
//#define EDGEFINDINGEXPL 6
//#define EDGEEXPL 5
//#define BOUNDEXPL 4
//#define CYCLEEXPL 3
//#define CLAUSEEXPL 2
//#define PATHEXPL 1
//#define NOEXPL 0


// arbitrary system to index both directions in a directed graph
#define OUT 0
#define IN 1


// the numeric gap between two values in numeric types
template <class T> class Gap {
public:
  static constexpr T epsilon() noexcept { return T(); }
};

template <> class Gap<short> {
public:
  static constexpr short epsilon() noexcept { return 1; }
};

template <> class Gap<int> {
public:
  static constexpr int epsilon() noexcept { return 1; }
};

template <> class Gap<long> {
public:
  static constexpr long int epsilon() noexcept { return 1; }
};

template <> class Gap<long long> {
public:
  static constexpr long long int epsilon() noexcept { return 1; }
};

template <> class Gap<unsigned short> {
public:
  static constexpr short epsilon() noexcept { return 1; }
};

template <> class Gap<unsigned int> {
public:
  static constexpr int epsilon() noexcept { return 1; }
};

template <> class Gap<unsigned long> {
public:
  static constexpr long int epsilon() noexcept { return 1; }
};

template <> class Gap<unsigned long long> {
public:
  static constexpr long long int epsilon() noexcept { return 1; }
};

template <> class Gap<float> {
public:
  static constexpr float epsilon() noexcept {
    return 1e-3; // std::numeric_limits<float>::epsilon();
  }
};

template <> class Gap<double> {
public:
  static constexpr double epsilon() noexcept {
    return 1e-6; // std::numeric_limits<double>::epsilon();
  }
};

template <> class Gap<long double> {
public:
  static constexpr double epsilon() noexcept {
    return 1e-9; // std::numeric_limits<long double>::epsilon();
  }
};

//
double cpu_time(void);

void seed(const unsigned long s);

unsigned long random(void);

} // namespace tempo

#endif // _TEMPO_SCHEDULING_HPP
