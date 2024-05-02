#ifndef __FORMAT_HH
#define __FORMAT_HH

#include <vector>

#include "Global.hpp"
#include "util/serialization.hpp"

template <typename T> class Resource : public std::vector<T> {

public:
  Resource() = default;

    Resource(std::vector<T> tasks, std::vector<T> demand, std::vector<std::vector<T>> transition, T capacity) :
            std::vector<T>(std::move(tasks)), demand(std::move(demand)), transition(std::move(transition)),
            capacity(capacity) {}

  T getTransitionTime(const size_t i, const size_t j) const {
    return (!transition.empty() ? transition[i][j] : 0);
  }
  T getDemand(const size_t i) const {
    return (!demand.empty() ? demand[i] : 1);
  }
  T getCapacity() const { return capacity; }
  //    task getTask(const size_t i) const {return task[i];}

  // private:

  //    std::vector<tempo::task> task;
  std::vector<T> demand;
  std::vector<std::vector<T>> transition;
  T capacity{1};
};

namespace nlohmann {
    template<typename T>
    struct adl_serializer<Resource<T>> {
        static void to_json(json &j, const Resource<T> &resource) {
            j["tasks"] = static_cast<const std::vector<T>&>(resource);
            j["demand"] = resource.demand;
            j["transition"] = resource.transition;
            j["capacity"] = resource.capacity;
        }

        static void from_json(const json &j, Resource<T> &resource) {
            resource = Resource<T>(j.at("tasks"), j.at("demand"), j.at("transition"), j.at("capacity"));
        }
    };

}

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
    
//    size_t nbTasks() { return durations.size(); }
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ProblemInstance, lowerBound, optimalSolution, durations, constraints, resources)


#endif
