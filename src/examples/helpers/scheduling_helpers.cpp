/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief contains utility that allows to easily load and solve scheduling problems
*/

#include <ranges>

#include "scheduling_helpers.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/parsing/osp.hpp"


auto loadSchedulingProblem(const tempo::Options &options) -> std::pair<SolverPtr, ProblemInstance> {
    using namespace tempo;
    using namespace std::views;
    auto solver = std::make_unique<Solver<>>(options);
    auto schedule{solver->newInterval(0, Constant::Infinity<int>, 0, 0, 0, Constant::Infinity<int>)};
    std::vector<DisjunctiveResource<>> resources;
    std::vector<Interval<>> tasks;

    if (options.input_format == "osp") {
        osp::parse(options.instance_file, *solver, schedule, tasks, resources);
    } else {
        std::cerr << "problem type " << options.input_format << " is not (yet) supported" << std::endl;
        std::exit(1);
    }


    for (const auto &consumingTasks: resources) {
        auto constraint = NoOverlap(schedule,
                                    consumingTasks | transform([&tasks](auto taskId) { return tasks.at(taskId); }),
                                    std::vector<std::vector<int>>{});
        solver->post(constraint);
    }

    SchedulingProblemHelper problem(std::move(tasks), std::move(resources), {}, schedule);
    auto trivialUb = 0;
    for (const auto &t : problem.tasks()) {
        trivialUb += t.minDuration(*solver);
    }

    trivialUb = std::min(trivialUb, options.ub);
    solver->set(schedule.end.before(trivialUb));
    return {std::move(solver), std::move(problem)};
}
