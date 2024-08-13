/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief contains utility that allows to easily load and solve scheduling problems
*/

#include <ranges>

#include "scheduling_helpers.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/jsp.hpp"


auto loadSchedulingProblem(const tempo::Options &options)
-> std::tuple<SolverPtr, ProblemInstance, std::optional<int>> {
    using Time = int;
    using namespace tempo;
    using namespace std::views;
    auto solver = std::make_unique<Solver<>>(options);
    auto schedule{solver->newInterval(0, Constant::Infinity<Time>, 0, 0, 0, Constant::Infinity<Time>)};
    std::vector<DisjunctiveResource<Time>> resources;
    std::vector<Interval<Time>> tasks;
    std::optional<Time> optSol;
    std::vector<DistanceConstraint<Time>> precedences;

    if (options.input_format == "osp") {
        optSol = osp::parse(options.instance_file, *solver, schedule, tasks, resources);
    } else if (options.input_format == "jsp") {
        jsp::parse(options.instance_file, *solver, schedule, tasks, resources, &precedences);
    } else {
        std::cerr << "problem type " << options.input_format << " is not (yet) supported" << std::endl;
        std::exit(1);
    }


    for (const auto &consumingTasks: resources) {
        auto constraint = NoOverlap(schedule,
                                    consumingTasks | transform([&tasks](auto taskId) { return tasks.at(taskId); }),
                                    std::vector<std::vector<Time>>{});
        solver->post(constraint);
    }

    SchedulingProblemHelper problem(std::move(tasks), std::move(resources), std::move(precedences), schedule);
    auto trivialUb = 0;
    for (const auto &t : problem.tasks()) {
        trivialUb += t.minDuration(*solver);
    }

    trivialUb = std::min(trivialUb, options.ub);
    solver->set(schedule.end.before(trivialUb));
    return {std::move(solver), std::move(problem), optSol};
}

void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch &branch) {
    for (auto [id, val] : branch) {
        auto lit = solver.boolean.getLiteral(val, id);
        solver.set(lit);
    }

    solver.propagate();
}
