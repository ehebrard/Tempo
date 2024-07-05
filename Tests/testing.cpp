/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <fstream>
#include "testing.hpp"

namespace tempo::testing {

    auto createTestProblem() -> ProblemInstance {
        using namespace tempo;
        std::vector<Interval<int>> tasks{{{0, 0}, {1, 0}, 2, 5},
                                         {{2, 0}, {2, 6}, 6, 6},
                                         {{3, 0}, {4, 0}, 2, 7},
                                         {{5, 0}, {6, 0}, 2, 3}};
        std::vector<Resource> resources{{2, {tasks[0], tasks[2]},           {2, 1}},
                                        {3, {tasks[1], tasks[2], tasks[3]}, {3, 1, 1}}};
        std::vector<DistanceConstraint<int>> precedences{{2, 1, 0}};
        return {std::move(tasks), std::move(resources), std::move(precedences)};
    }

    auto createTasks(std::initializer_list<std::pair<int, int>> durations) -> std::vector<Interval<int>> {
        std::vector<Interval<int>> ret;
        ret.reserve(durations.size());
        unsigned idx = 0;
        for (auto [min, max] : durations) {
            if (min == max) {
                ret.emplace_back(TemporalVar{idx, 0}, TemporalVar{idx, min}, min, min);
            } else {
                ret.emplace_back(TemporalVar{idx, 0}, TemporalVar{idx + 1, 0}, min, max);
                ++idx;
            }

            ++idx;
        }

        return ret;
    }

    auto createDummyTasks(unsigned int numberOfTasks) -> std::vector<Interval<int>> {
        std::vector<Interval<int>> ret;
        ret.reserve(numberOfTasks);
        for (unsigned i = 0; i < numberOfTasks; ++i) {
            ret.emplace_back(TemporalVar{i, 0}, TemporalVar{i, 0}, 0, 0);
        }

        return ret;
    }

    Resource::Resource(int capacity, std::vector<Interval<int>> tasks, std::vector<int> demands)
            : std::vector<Interval<int>>(std::move(tasks)), demands(std::move(demands)), capacity(capacity) {}

    int Resource::getDemand(unsigned int taskId) const {
        return demands.at(taskId);
    }

    int Resource::resourceCapacity() const {
        return capacity;
    }
}

