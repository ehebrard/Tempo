
#ifndef _TEMPO_OBJECTIVE_HPP
#define _TEMPO_OBJECTIVE_HPP


#include "util/parsing/format.hpp"
//#include "Scheduler.hpp"

namespace tempo {

template <typename T> class Scheduler;

template <typename T> class Objective {
public:
  Objective() {}
  ~Objective() = default;

  T gap() { return ub - lb; }
  void closeGap() { lb = ub; }
  T lowerBound() const { return lb; }
  T upperBound() const { return ub; }

  std::ostream &display(std::ostream &os) const {
    os << "[" << std::left << std::setw(5) << std::setfill('.') << lowerBound()
       << std::setfill(' ');
    auto ub{upperBound()};
    if (ub < INFTY)
      os << std::right << std::setw(6) << std::setfill('.') << ub
         << std::setfill(' ');
    else
      os << ".infty";
    os << "]";
    return os;
  }

protected:
  T lb;
  T ub;
};

template <typename T> class Makespan : public Objective<T> {
public:
  Makespan(Scheduler<T> &s) : Objective<T>(), schedule(s) {}
  ~Makespan() = default;

  void updateLowerBound() { Objective<T>::lb = schedule.lower(HORIZON); }

  void updateUpperBound() { Objective<T>::ub = schedule.lower(HORIZON); }
    
    void applyUpperBound() { schedule.newMaximumLag(ORIGIN, HORIZON, Objective<T>::ub - Gap<T>::epsilon()); }

private:
  Scheduler<T> &schedule;
};


template <typename T> class PathLength : public Objective<T> {
public:
    PathLength(Scheduler<T> &s, Resource<T>& r) : Objective<T>(), schedule(s), stops(r) { Objective<T>::lb = 0; }
  ~PathLength() = default;

    void updateLowerBound() {
        if(Objective<T>::lb == 0) {
            for(auto x : stops) {
                Objective<T>::lb += schedule.minDuration(x);
            }
            for(auto& transitions : stops.transition) {
                Objective<T>::lb += *(std::min_element(transitions.begin(), transitions.end()));
            }
        }
    }

    void updateUpperBound() {
        std::sort(stops.begin(), stops.end(), [&](const task a, const task b) {return schedule.lower(START(a)) < schedule.lower(START(b));} );
        
        Objective<T>::ub = schedule.minDuration(stops[0]);
        for(size_t i{1}; i<stops.size(); ++i) {
            Objective<T>::ub += stops.transition[stops[i-1]][stops[i]];
            Objective<T>::ub += schedule.minDuration(stops[i]);
        }
//        Objective<T>::ub = schedule.lower(HORIZON);
    }
    
    void applyUpperBound() { std::cout << "TODO!\n"; }

private:
  Scheduler<T> &schedule;
    Resource<T> &stops;
};

template <typename T> class NoObj : public Objective<T> {
public:
    NoObj() : Objective<T>() {Objective<T>::lb = 0;}
  ~NoObj() = default;

  void updateLowerBound() {  }

  void updateUpperBound() { Objective<T>::ub = 1; }
};

template <typename T> class No {
public:
  static NoObj<T> Obj;
};

template <typename T> NoObj<T> No<T>::Obj = NoObj<T>();

} // namespace tempo

#endif
