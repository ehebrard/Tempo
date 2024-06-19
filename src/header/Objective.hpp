
#ifndef _TEMPO_OBJECTIVE_HPP
#define _TEMPO_OBJECTIVE_HPP

//#include "util/parsing/format.hpp"
//#include "constraints/Cardinality.hpp"

namespace tempo {

template <typename T> class Scheduler;

//template <typename T> class Solver;

template <typename T> class Job;



template <typename T> class MakespanObjective {
public:
    MakespanObjective(Job<T> &j, Solver<T>& s) : job(j), solver(s) {}
    MakespanObjective(Job<T> &j, Solver<T>& s, const T u) : job(j), solver(s) { setPrimal(u); }
  ~MakespanObjective() = default;

  T value() { return job.getEarliestEnd(solver); }

  T gap() { return p_b - d_b; }
  T dualBound() const { return d_b; }
  T primalBound() const { return p_b; }

  void setDual(const T v) { d_b = v; }
  void initDual() { d_b = job.getEarliestEnd(solver); }

  void setPrimal(const T v) {
    p_b = v;
    if (gap()) {
      apply(p_b - Gap<T>::epsilon());
    }
  }

  void apply(const T target) {
      solver.set(job.end.before(target));
  }

  std::ostream &display(std::ostream &os) const {
    os << "[" << std::left << std::setw(5) << std::setfill('.') << dualBound()
       << std::setfill(' ');
    auto pb{primalBound()};
    if (pb < Constant::Infinity<T>)
      os << std::right << std::setw(6) << std::setfill('.') << pb
         << std::setfill(' ');
    else
      os << ".infty";
    os << "]";
    return os;
  }

private:
  Job<T> &job;
    Solver<T> &solver;
  T d_b{0};
    T p_b{Constant::Infinity<T>};
};

//template <typename T> class MaximumCardinality {
//public:
//  template <typename Iter>
//  MaximumCardinality(Scheduler<T> &s, Iter beg_lit, Iter end_lit)
//      : schedule(s) {
//    for (auto l{beg_lit}; l != end_lit; ++l)
//      literals.push_back(*l);
//    card = new CardinalityGeq<T>(schedule, beg_lit, end_lit, 0);
//    schedule.post(card);
//  }
//  ~MaximumCardinality() = default;
//
//  T gap() { return p_b - d_b; }
//  //  void closeGap() { d_b = p_b; }
//  T dualBound() const { return d_b; }
//  T primalBound() const { return p_b; }
//
//  size_t value() {
//    size_t c{0};
//    for (auto l : literals) {
//      c += schedule.satisfied(l);
//    }
//    return c;
//  }
//
//  void setDual(const T v) { d_b = v; }
//
//  void setPrimal(const T v) {
//    p_b = v;
//    if (gap()) {
//      apply(p_b + Gap<T>::epsilon());
//    }
//  }
//
//  void apply(const int target) { card->setBound(target); }
//
//private:
//  Scheduler<T> &schedule;
//  std::vector<lit> literals;
//  CardinalityGeq<T> *card;
//  T d_b;
//  T p_b;
//};

template <typename T> class NoObj {
public:
  NoObj() = default;
  ~NoObj() = default;

  T gap() { return 1; }
  T dualBound() const { return 0; }

  void setDual(const T) {}

  T value() { return 0; }

  std::ostream &display(std::ostream &os) const {
    os << " no solution ";
    return os;
  }
};

template <typename T> class No {
public:
  static NoObj<T> Obj;
};

template <typename T> NoObj<T> No<T>::Obj = NoObj<T>();

} // namespace tempo

#endif
