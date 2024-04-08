#ifndef __FORMAT_HH
#define __FORMAT_HH

#include <vector>

#include "Global.hpp"

template <typename T> class Resource : public std::vector<T> {

public:
  Resource() = default;
  ~Resource() = default;

  T getTransitionTime(const size_t i, const size_t j) const {
    return (transition.size() > 0 ? transition[i][j] : 0);
  }
  T getDemand(const size_t i) const {
    return (demand.sizee() > 0 ? demand[i] : 1);
  }
  T getCapacity() const { return capacity; }
  //    task getTask(const size_t i) const {return task[i];}

  // private:

  //    std::vector<tempo::task> task;
  std::vector<T> demand;
  std::vector<std::vector<T>> transition;
  T capacity{1};
};

// template<typename T>
// class Resource {
//
// public:
//     Resource() = default;
//     ~Resource() = default;
//
//     virtual T getTransitionTime(const size_t i, const size_t j) const = 0;
//     virtual T getDemand(const size_t i) const = 0;
//     virtual T getCapacity() const = 0;
//     task getTask(const size_t i) const = 0;
//
// private:
//
//     std::vector<tempo::task> tasks;
//
// };

// template<typename T>
// class DisjunctiveResource : public Resource<T> {
//
// public:
//     Resource() = default;
//     ~Resource() = default;
//
//     virtual T getTransitionTime(const size_t i, const size_t j) const = 0;
//     virtual T getDemand(const size_t i) const = 0;
//     virtual T getCapacity() const = 0;
//     task getTask(const size_t i) const = 0;
//
// private:
//
//     std::vector<tempo::task> tasks;
//
// };

struct ProblemInstance {
    int lowerBound;
    int optimalSolution;
    std::vector<int> durations;
    std::vector<std::tuple<tempo::event, tempo::event, int>> constraints;
    std::vector<Resource<int>> resources;
    //    std::vector<std::vector<tempo::task>> resources;
    //    std::vector<std::vector<int>> transition;
};


#endif
