/**
 * @author Tim Luchterhand
 * @date 27.06.23.
 */

#include <fstream>
#include "testing.hpp"

namespace tempo::testing {

    auto createTestProblem() -> std::pair<ProblemInstance, DummyScheduler> {
        using namespace tempo;
        auto [tasks, scheduler] = createTasks({TaskSpec{.minDur = 2, .maxDur = 5},
                                               TaskSpec{.minDur = 6, .maxDur = 6},
                                               TaskSpec{.minDur = 2, .maxDur = 7},
                                               TaskSpec{.minDur = 2, .maxDur = 3},
                                               TaskSpec{.minDur = 0, .maxDur = Constant::Infinity<int>}});
        auto schedule = tasks.back();
        tasks.pop_back();
        std::vector<Resource> resources{{2, {tasks[0], tasks[2]},           {2, 1}},
                                        {3, {tasks[1], tasks[2], tasks[3]}, {3, 1, 1}}};
        std::vector<DistanceConstraint<int>> precedences{{tasks[1].start.id(), tasks[0].end.id(), 0}};
        return {ProblemInstance(std::move(tasks), std::move(resources), std::move(precedences), schedule),
                std::move(scheduler)} ;
    }

    auto createExtendedTestProblem() -> std::pair<ProblemInstance, DummyScheduler> {
        auto [problem, scheduler] = createTestProblem();
        auto &resources = const_cast<std::remove_cvref_t<decltype(problem.resources())>&>(problem.resources());
        resources.emplace_back(2, std::vector{problem.tasks().at(1), problem.tasks().at(3)}, std::vector{2, 1});
        return {std::move(problem), std::move(scheduler)};
    }


    auto createTasks(const std::vector<TaskSpec> &specs) -> std::pair<std::vector<Interval<int>>, DummyScheduler> {
        std::vector<Interval<int>> ret;
        std::vector<int> upper(specs.size() * VarTaskMapping::NumTemporalVarPerTask, Constant::Infinity<int>);
        std::vector<int> lower(specs.size() * VarTaskMapping::NumTemporalVarPerTask, 0);
        unsigned idx = 0;
        for (const auto &spec : specs) {
            if (spec.minDur == spec.maxDur) {
                ret.emplace_back(NumericVar{idx, 0}, NumericVar{idx, spec.minDur}, NumericVar{idx + 1, 0});
            } else {
                ret.emplace_back(NumericVar{idx, 0}, NumericVar{idx + 1, 0}, NumericVar{idx + 2, 0});
                ++idx;
            }

            lower.at(ret.back().duration.id()) = spec.minDur;
            upper.at(ret.back().duration.id()) = spec.maxDur;
            lower.at(ret.back().start.id()) = spec.earliestStart;
            upper.at(ret.back().end.id()) = spec.latestDeadline;
            idx += 2;
        }

        return {std::move(ret), DummyScheduler(std::move(upper), std::move(lower))};
    }

    auto createDummyTasks(unsigned int numberOfTasks) -> std::vector<Interval<int>> {
        std::vector<Interval<int>> ret;
        ret.reserve(numberOfTasks);
        for (unsigned i = 0; i < numberOfTasks; ++i) {
            ret.emplace_back(NumericVar{i, 0}, NumericVar{i, 0}, NumericVar{i, 0});
        }

        return ret;
    }

    auto createRandomProblem(std::size_t numTasks, std::size_t numResources,
                             double precedenceChance) -> std::tuple<ProblemInstance, DummyScheduler, Matrix<int>> {
        assert(precedenceChance >= 0 and precedenceChance <= 1);
        using namespace tempo;
        std::vector<DistanceConstraint<int>> precedences;
        Matrix<int> eventNetwork(2 * numTasks + 2, 2 * numTasks + 2);
        int upperBound = 0;
        std::vector<TaskSpec> taskSpecs;
        taskSpecs.reserve(numTasks + 1);
        for (int t = 0; t < static_cast<int>(numTasks); ++t) {
            taskSpecs.emplace_back(random_int(2, 6), random_int(7, 10), random_int(-7, 0), random_int(-7, 0));
            upperBound += taskSpecs.back().maxDur;
        }

        taskSpecs.emplace_back(upperBound / 2, upperBound, 0, upperBound);
        auto [tasks, scheduler] = createTasks(taskSpecs);
        auto sched = tasks.back();
        tasks.pop_back();
        Matrix<int> taskDistances(numTasks, numTasks, upperBound);
        for (unsigned from = 0; from < numTasks; ++from) {
            for(unsigned to = 0; to < numTasks; ++to) {
                if (from == to) {
                    continue;
                }

                taskDistances.at(from, to) = random_int(-20 * static_cast<int>(precedenceChance),
                                                        20 * static_cast<int>(1 - precedenceChance));
                if (taskDistances.at(from, to) <= 0) {
                    precedences.emplace_back(tasks.at(from).start.id(), tasks.at(to).end.id(),
                                             eventNetwork.at(from, to));
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

        return {ProblemInstance(std::move(tasks), std::move(resources), std::move(precedences), sched),
                std::move(scheduler), std::move(taskDistances)};

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

