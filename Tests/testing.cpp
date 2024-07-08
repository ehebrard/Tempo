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
        return {std::move(tasks), std::move(resources), std::move(precedences),
                {{7, 0}, {8, 0}, 0, Constant::Infinity<int>}};
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

    auto createRandomProblem(std::size_t numTasks, std::size_t numResources,
                             double precedenceChance) -> std::tuple<ProblemInstance, DummyScheduler, Matrix<int>> {
        assert(precedenceChance >= 0 and precedenceChance <= 1);
        using namespace tempo;
        std::vector<Interval<int>> tasks;
        std::vector<DistanceConstraint<int>> precedences;
        tasks.reserve(numTasks);
        Matrix<int> eventNetwork(2 * numTasks + 2, 2 * numTasks + 2);
        int upperBound = 0;
        std::vector<int> upper(numTasks * 2 + 2, -1);
        std::vector<int> lower(numTasks * 2 + 2, -1);
        for (int t = 0; t < static_cast<int>(numTasks); ++t) {
            tasks.emplace_back(TemporalVar(t, 0), TemporalVar(t + 1, 0), random_int(2, 6), random_int(7, 10));
            lower.at(tasks.back().start) = random_int(-7, 0);
            upper.at(tasks.back().end) = random_int(-7, 0);
            upperBound += upper.at(tasks.back().end);
        }

        const Interval<int> sched({2 * static_cast<var_t>(numTasks), 0}, {2 * static_cast<var_t>(numTasks) + 1, 0}, 0,
                                  Constant::Infinity<int>);
        upper.at(sched.end) = upperBound;
        lower.at(sched.end) = upperBound / 2;
        Matrix<int> taskDistances(numTasks, numTasks, upperBound);

        for (unsigned from = 0; from < numTasks; ++from) {
            for(unsigned to = 0; to < numTasks; ++to) {
                if (from == to) {
                    continue;
                }

                taskDistances.at(from, to) = random_int(-20 * static_cast<int>(precedenceChance),
                                                        20 * static_cast<int>(1 - precedenceChance));
                if (taskDistances.at(from, to) <= 0) {
                    precedences.emplace_back(tasks.at(from).start, tasks.at(to).end, eventNetwork.at(from, to));
                }
            }
        }

        std::vector<Resource> resources;
        resources.reserve(numResources);
        for (auto r = 0ul; r < numResources; ++r) {
            decltype(tasks) consuming;
            std::vector<int> demands;
            const auto capacity = random_int(1, 5);
            for (auto t = 0ul; t < numTasks; ++t) {
                auto demand = random_int(0, capacity);
                if (demand > 0) {
                    consuming.emplace_back(tasks.at(t));
                    demands.emplace_back(demand);
                }
            }

            resources.emplace_back(capacity, std::move(consuming), std::move(demands));
        }

        return {{std::move(tasks), std::move(resources), std::move(precedences), sched},
                DummyScheduler(std::move(upper), std::move(lower)), std::move(taskDistances)};

    }

    Resource::Resource(int capacity, std::vector<Interval<int>> tasks, std::vector<int> demands)
            : std::vector<Interval<int>>(std::move(tasks)), demands(std::move(demands)), capacity(capacity) {}

    int Resource::getDemand(unsigned int taskId) const {
        return demands.at(taskId);
    }

    int Resource::resourceCapacity() const {
        return capacity;
    }

    BoundProvider::BoundProvider(std::vector<int> upper, std::vector<int> lower) : u(std::move(upper)), l(std::move(lower)) {}

    int BoundProvider::upper(tempo::var_t var) const { return u.at(var); }

    int BoundProvider::lower(tempo::var_t var) const { return l.at(var); }
}

