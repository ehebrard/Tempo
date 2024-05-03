
#ifndef _TEMPO_OBJECTIVE_HPP
#define _TEMPO_OBJECTIVE_HPP

//#include "util/parsing/format.hpp"
#include "constraints/Cardinality.hpp"

namespace tempo {

template <typename T> class Scheduler;

// template <typename T> class Objective {
// public:
////  Objective() {}
////  Objective(const T u) : p_b(u) {}
////  ~Objective() = default;
////
////  T gap() { return p_b - d_b; }
////  //  void closeGap() { d_b = p_b; }
////  T dualBound() const { return d_b; }
////  T primalBound() const { return p_b; }
//////  void setDual(const T v) { d_b = v; }
//////  void setPrimal(const T v) { p_b = v; }
//
//  std::ostream &display(std::ostream &os) const {
//    os << "[" << std::left << std::setw(5) << std::setfill('.') << dualBound()
//       << std::setfill(' ');
//    auto pb{primalBound()};
//    if (pb < INFTY)
//      os << std::right << std::setw(6) << std::setfill('.') << pb
//         << std::setfill(' ');
//    else
//      os << ".infty";
//    os << "]";
//    return os;
//  }
//
////protected:
////  T d_b;
////  T p_b;
//};

template <typename T> class Makespan {
public:
  Makespan(Scheduler<T> &s) : schedule(s) {}
  Makespan(Scheduler<T> &s, const T u) : schedule(s) { setPrimal(u); }
  ~Makespan() = default;

  //  void updatedualBound() { Objective<T>::d_b = schedule.lower(HORIZON); }
  //
  //  void updateprimalBound() { Objective<T>::p_b = schedule.lower(HORIZON); }
  T value() { return schedule.lower(HORIZON); }

  T gap() { return p_b - d_b; }
  //  void closeGap() { d_b = p_b; }
  T dualBound() const { return d_b; }
  T primalBound() const { return p_b; }

  void setDual(const T v) { d_b = v; }
  void initDual() { d_b = schedule.lower(HORIZON); }

  void setPrimal(const T v) {
    p_b = v;
    if (gap()) {
      apply(p_b - Gap<T>::epsilon());
    }
  }

  void apply(const T target) {
    schedule.newMaximumLag(ORIGIN, HORIZON, target);
  }

  std::ostream &display(std::ostream &os) const {
    os << "[" << std::left << std::setw(5) << std::setfill('.') << dualBound()
       << std::setfill(' ');
    auto pb{primalBound()};
    if (pb < INFTY)
      os << std::right << std::setw(6) << std::setfill('.') << pb
         << std::setfill(' ');
    else
      os << ".infty";
    os << "]";
    return os;
  }

private:
  Scheduler<T> &schedule;
  T d_b{0};
  T p_b;
};

template <typename T> class MaximumCardinality {
public:
  template <typename Iter>
  MaximumCardinality(Scheduler<T> &s, Iter beg_lit, Iter end_lit)
      : schedule(s) {
    for (auto l{beg_lit}; l != end_lit; ++l)
      literals.push_back(*l);
    card = new CardinalityGeq<T>(schedule, beg_lit, end_lit, 0);
    schedule.post(card);
  }
  ~MaximumCardinality() = default;

  T gap() { return p_b - d_b; }
  //  void closeGap() { d_b = p_b; }
  T dualBound() const { return d_b; }
  T primalBound() const { return p_b; }

  size_t value() {
    size_t c{0};
    for (auto l : literals) {
      c += schedule.satisfied(l);
    }
    return c;
  }

  void setDual(const T v) { d_b = v; }

  void setPrimal(const T v) {
    p_b = v;
    if (gap()) {
      apply(p_b + Gap<T>::epsilon());
    }
  }

  void apply(const int target) { card->setBound(target); }

private:
  Scheduler<T> &schedule;
  std::vector<lit> literals;
  CardinalityGeq<T> *card;
  T d_b;
  T p_b;
};

// template <typename T> class PathLength : public Objective<T> {
// public:
//     PathLength(Scheduler<T> &s, Resource<T>& r) : Objective<T>(),
//     schedule(s), stops(r) { Objective<T>::        = 0; }
//   ~PathLength() = default;
//
//   //    void updatedualBound() {
//   //        if(Objective<T>::lb == 0) {
//   //            for(auto x : stops) {
//   //                Objective<T>::d_b += schedule.minDuration(x);
//   //            }
//   //            for(auto& transitions : stops.transition) {
//   //                Objective<T>::d_b +=
//   *(std::min_element(transitions.begin(),
//   //                transitions.end()));
//   //            }
//   //        }
//   //    }
//   //
//   //    void updateprimalBound() {
//   //        std::sort(stops.begin(), stops.end(), [&](const task a, const
//   task
//   //        b) {return schedule.lower(START(a)) < schedule.lower(START(b));}
//   );
//   //
//   //        Objective<T>::ub = schedule.minDuration(stops[0]);
//   //        for(size_t i{1}; i<stops.size(); ++i) {
//   //            Objective<T>::ub += stops.transition[stops[i-1]][stops[i]];
//   //            Objective<T>::ub += schedule.minDuration(stops[i]);
//   //        }
//   ////        Objective<T>::ub = schedule.lower(HORIZON);
//   //    }
//   T value() {
//
//     std::sort(stops.begin(), stops.end(), [&](const task a, const task b) {
//       return schedule.lower(START(a)) < schedule.lower(START(b));
//     });
//
//     T val{schedule.minDuration(stops[0])};
//     for (size_t i{1}; i < stops.size(); ++i) {
//       val += stops.transition[stops[i - 1]][stops[i]];
//       val += schedule.minDuration(stops[i]);
//     }
//     return val;
//   }
//
//   void apply(const T) { std::cout << "TODO!\n"; }
//
// private:
//   Scheduler<T> &schedule;
//     Resource<T> &stops;
// };

template <typename T> class NoObj {
public:
  NoObj() = default;
  ~NoObj() = default;

  //  void updatedualBound() {  }
  //
  //  void updateprimalBound() { Objective<T>::ub = 1; }

  T gap() { return 1; }
  T dualBound() const { return 0; }
  //    T primalBound() const { return 1; }

  void setDual(const T) {}

  T value() { return 0; }

  std::ostream &display(std::ostream &os) const {
    os << " no solution ";
    return os;
  }

  //    void updateprimalBound() { Objective<T>::ub = 1; }
};

template <typename T> class No {
public:
  static NoObj<T> Obj;
};

template <typename T> NoObj<T> No<T>::Obj = NoObj<T>();

} // namespace tempo

#endif
