#ifndef _TEMPO_GLOBAL_HPP
#define _TEMPO_GLOBAL_HPP

#include <string>
#include <limits>
#include <cstdint>
#include <sys/resource.h>


namespace tempo {

using var_t = std::uint32_t;
using info_t = std::uint32_t;

//#define DBG_BOUND (num_fails >= 20299)
//#define DBG_CBOUND (solver.num_fails >= 20299)
//#define DBG_BBOUND (sched.num_fails >= 20299) //(sched.num_fails >= 236)
//#define DBG_TRACE 33                             // 183 //1+2+4+32+128
#define SEARCH 1
#define DOMAINS 2
#define BRANCH 4
#define CLAUSES 8
#define PROPAGATION 16
#define LEARNING 32
#define QUEUE 64
#define UNITPROPAGATION 128

//#define DBG_MINIMIZATION
//#define DBG_EDGEFINDING (m_schedule.num_cons_propagations >= 55553)
//#define DBG_EXPLEF true      //(m_schedule.num_fails > 236)
//#define DBG_THETA (m_schedule.num_fails >= 481)
//#define DBG_BELLMAN true //(sched.num_choicepoints >= 1045)
//#define DBG_BELLMAN_EXPL (sched.num_choicepoints >= 4172)
//#define DEBUG_HEURISTICS
//#define DBG_UP
//#define DBG_CL 20299
//#define
//#define DBG_CLPLUS 1030
//#define DBG_TRANSITIVITY true //(m_schedule.num_choicepoints >= 4064)
//#define DBG_EXPL_TRANS true
//#define DBG_SOL

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

// using index_t = size_t;
using index_t = std::uint32_t;

using event = int;
#define NOEVENT -1
#define ORIGIN 0
#define HORIZON 1

using task = int;

using var = int;
#define NoVar -1

using genlit = int;
using lit = int;
#define NoLit -1

using boolean_state = int;
#define True 1
#define False 0
#define Unknown -1

using hint = int;
#define NoHint -1

using lit_type = int;
#define BOUND_LIT 0
#define EDGE_LIT 1

#define CARDEXPL 7
#define TRANSITIVITYEXPL 7
#define EDGEFINDINGEXPL 6
#define EDGEEXPL 5
#define BOUNDEXPL 4
#define CYCLEEXPL 3
#define CLAUSEEXPL 2
#define PATHEXPL 1
#define NOEXPL 0


//#define POSITIVE 1
//#define NEGATIVE 0

#define INFTY std::numeric_limits<int32_t>::max()

#define OUT 0
#define IN 1

#define LOWER 0
#define UPPER 1


task TASK(event x);
event START(task x);
event END(task x);

lit LIT(const var x, boolean_state v);
lit POS(const var x);
lit NEG(const var x);
lit NOT(const lit l);
genlit GNOT(const genlit l);

var VAR(const lit l);

lit LOWERBOUND(const event x);
lit UPPERBOUND(const event x);

event EVENT(const lit l);

boolean_state SIGN(const lit l);

lit_type LTYPE(const genlit h);
genlit BOUND(const lit l);
genlit EDGE(const lit l);
lit FROM_GEN(const genlit h);

int TODIMACS(const lit l);


template<typename T>
bool finite(const T x) {
    return x < INFTY/4;
}

char etype(const event x);
std::string prettyEvent(const event e);
std::string prettyEventLit(const lit el);
//std::string prettyLiteral(const genlit el);


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

double cpu_time(void);

void seed(const unsigned long s);

unsigned long random(void);

} // namespace tempo

#endif // _TEMPO_SCHEDULING_HPP
