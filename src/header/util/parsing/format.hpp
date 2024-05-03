#ifndef __FORMAT_HH
#define __FORMAT_HH

#include <vector>

#include "Global.hpp"

template <typename T> class Resource : public std::vector<int> {

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


struct ProblemInstance {
    
    ProblemInstance() = default;
    ~ProblemInstance() = default;
    
    int lowerBound;
    int optimalSolution;
    std::vector<int> durations;
    std::vector<std::tuple<tempo::event, tempo::event, int>> constraints;
    std::vector<Resource<int>> resources;
    //    std::vector<std::vector<tempo::task>> resources;
    //    std::vector<std::vector<int>> transition;
    
//    size_t nbTasks() { return durations.size(); }
    
    
};


#endif
