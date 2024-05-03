/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <fstream>
#include "testing.hpp"

namespace tempo::testing {
    auto createRandomProblem(std::size_t numTasks, std::size_t numResources,
                             double precedenceChance) -> std::pair<ProblemInstance, Matrix<int>> {
        assert(precedenceChance >= 0 and precedenceChance <= 1);
        std::vector<int> durations;
        std::vector<std::tuple<event, event, int>> precedences;
        durations.reserve(numTasks);
        Matrix<int> eventNetwork(2 * numTasks + 2, 2 * numTasks + 2);
        int upperBound = 0;
        for (task t = 0; t < static_cast<task>(numTasks); ++t) {
            TaskSpec spec{random_int(2, 6), random_int(6, 10), random_int(-7, 0), random_int(-7, 0)};
            durations.emplace_back(spec.maxDur);
            upperBound += spec.maxDur;
            setTaskDurations(t, spec.minDur, spec.maxDur, spec.release, spec.deadline, eventNetwork);
        }

        setUpperBound(upperBound, eventNetwork);
        setLowerBound(upperBound / 2, eventNetwork);

        for (event from = START(0); from < END(static_cast<int>(numTasks) - 1); ++from) {
            for(event to = START(0); to < END(static_cast<int>(numTasks) - 1); ++to) {
                if (from == to) {
                    continue;
                }

                eventNetwork.at(from, to) = random_int<int>(-20 * precedenceChance, 20 * (1 - precedenceChance));
                if (eventNetwork.at(from, to) <= 0) {
                    precedences.emplace_back(from, to, eventNetwork.at(from, to));
                }
            }
        }

        std::vector<Resource<int>> resources;
        resources.reserve(numResources);
        for (std::size_t r = 0; r < numResources; ++r) {
            const auto capacity = random_int(0, 5);
            std::vector<int> demands;
            std::vector<int> tasksIds;
            for (std::size_t t = 0; t < numTasks; ++t) {
                auto demand = random_int(0, capacity);
                if (demand > 0) {
                    demands.emplace_back(demand);
                    tasksIds.emplace_back(t);
                }
            }

            resources.emplace_back(std::move(tasksIds), std::move(demands), std::vector<std::vector<int>>{}, capacity);
        }

        return {ProblemInstance{.lowerBound = upperBound / 2, .optimalSolution = upperBound / 2,
                                .durations = std::move(durations), .constraints = std::move(precedences),
                                .resources = std::move(resources)}, std::move(eventNetwork)};
    }
}

