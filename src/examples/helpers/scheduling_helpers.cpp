/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief contains utility that allows to easily load and solve scheduling problems
*/

#include <ranges>
#include <Iterators.hpp>

#include "scheduling_helpers.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/parsing/osp.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/jstl.hpp"
#include "util/parsing/rcpsp.hpp"


static auto parseDisjunctive(const tempo::Options &options, tempo::Solver<tempo::DefaultTime> &solver,
                             tempo::Interval<tempo::DefaultTime> schedule)
-> std::optional<std::pair<ProblemInstance, std::optional<tempo::DefaultTime>>> {
    std::optional<tempo::DefaultTime> optSol;
    std::vector<tempo::Interval<tempo::DefaultTime>> tasks;
    std::vector<DisjunctiveResource> resources;
    std::vector<tempo::DistanceConstraint<tempo::DefaultTime>> precedences;
    if (options.input_format == "osp") {
        optSol = osp::parse(options.instance_file, solver, schedule, tasks, resources);
    } else if (options.input_format == "jsp") {
        jsp::parse(options.instance_file, solver, schedule, tasks, resources, &precedences);
    } else if (options.input_format == "jstl") {
        jstl::parse(options.instance_file, solver, schedule, tasks, resources, &precedences);
    } else {
        return {};
    }

    auto tview = resources | std::views::transform([](auto &&r) { return Resource(std::move(r)); });
    std::vector<Resource> res(tview.begin(), tview.end());
    return {{ProblemInstance(std::move(tasks), std::move(res), std::move(precedences), schedule), optSol}};
}

static auto parseCumulative(const tempo::Options &options, tempo::Solver<tempo::DefaultTime> &solver,
                            tempo::Interval<tempo::DefaultTime> schedule)
-> std::optional<std::pair<ProblemInstance, std::optional<tempo::DefaultTime>>> {
    std::vector<tempo::Interval<tempo::DefaultTime>> tasks;
    std::vector<std::vector<std::size_t>> resources; // elems are the resource requirements of a task
    std::vector<std::vector<int>> demands; // elems are the resource demands of a task
    std::vector<std::pair<int, int>> p;
    std::vector<ResourceUnit> capacities; // resource capacities
    std::vector<tempo::DistanceConstraint<tempo::DefaultTime>> precedences;
    if (options.input_format == "rcpsp") {
        rcpsp::parse(options.instance_file, solver, schedule, tasks, resources, demands, capacities, p, &precedences);
    } else {
        return {};
    }
    const auto numResources = capacities.size();
    std::vector<std::vector<unsigned>> tasksByResource(numResources);
    std::vector<std::vector<ResourceUnit>> demandsByResource(numResources);
    for (auto [taskId, reqResources, resDemands]: iterators::const_zip_enumerate(resources, demands)) {
        for (auto [resId, resDem]: iterators::zip(reqResources, resDemands)) {
            tasksByResource.at(resId).emplace_back(taskId);
            demandsByResource.at(resId).emplace_back(resDem);
        }
    }

    std::vector<Resource> cResources;
    cResources.reserve(numResources);
    for (auto [taskList, demandsList, capacity]: iterators::zip(tasksByResource, demandsByResource, capacities)) {
        cResources.emplace_back(std::in_place_type<CumulativeResource<ResourceUnit>>, std::move(taskList),
                                std::move(demandsList), capacity);
    }

    return {{ProblemInstance(std::move(tasks), std::move(cResources), std::move(precedences), schedule), {}}};
}

static void postDisjunctive(tempo::Solver<tempo::DefaultTime> &solver, const ProblemInstance &problem) {
    const auto &tasks = problem.tasks();
    for (const auto &consumingTasks: problem.resources()) {
        auto constraint = tempo::NoOverlap(problem.schedule(),
                                           consumingTasks.tasks() |
                                           std::views::transform([&tasks](auto taskId) { return tasks.at(taskId); }),
                                           std::vector<std::vector<tempo::DefaultTime>>{});
        solver.post(constraint);
    }
}

void postCumulative(tempo::Solver<tempo::DefaultTime> &solver, const ProblemInstance &problem) {
    using namespace std::views;
    using namespace tempo;
    for (const auto &resource: problem.resources()) {
        auto capacity = solver.newConstant(resource.resourceCapacity());
        auto tView = resource.tasks() | transform([&problem](auto taskId) { return problem.tasks().at(taskId); });
        auto demView = iota(0u, static_cast<unsigned>(resource.tasks().size())) |
                       transform([&resource, &solver](auto idx) {
                           return solver.newConstant(resource.getDemand(idx));
                       });
        std::vector<Interval<tempo::DefaultTime>> tasks;
        tasks.reserve(tView.size());
        std::vector<NumericVar<ResourceUnit>> demands;
        demands.reserve(demView.size());
        for (auto [t, d]: iterators::zip(tView, demView)) {
            tasks.emplace_back(t);
            demands.emplace_back(d);
        }

        solver.post(Cumulative(problem.schedule(), capacity, tasks, demands));
    }
}

auto loadSchedulingProblem(const tempo::Options &options)
-> std::tuple<SolverPtr, ProblemInstance, std::optional<tempo::DefaultTime>, unsigned> {
    using namespace tempo;
    using namespace std::views;
    auto solver = std::make_unique<Solver<>>(options);
    const auto schedule{solver->newInterval(0, Constant::Infinity<tempo::DefaultTime>, 0, 0, 0,
                                            Constant::Infinity<tempo::DefaultTime>)};
    auto res = parseDisjunctive(options, *solver, schedule);
    if (res.has_value()) {
        auto &[problem, optSol] = *res;
        postDisjunctive(*solver, problem);
    } else if ((res = parseCumulative(options, *solver, schedule)).has_value()) {
        auto [problem, _] = *res;
        postCumulative(*solver, problem);
    } else {
        std::cerr << "problem type " << options.input_format << " is not (yet) supported" << std::endl;
        std::exit(1);
    }

    auto &[problem, optSol] = *res;
    tempo::DefaultTime trivialUb = 0;
    for (const auto &t: problem.tasks()) {
        auto maxDur = t.maxDuration(*solver);
        if (maxDur == Constant::Infinity<tempo::DefaultTime>) {
            trivialUb = maxDur;
            break;
        }

        trivialUb += maxDur;
    }

    trivialUb = std::min(trivialUb, options.ub);
    solver->set(schedule.end.before(trivialUb));
    auto numTasks = static_cast<unsigned>(problem.tasks().size());
    return {std::move(solver), std::move(problem), optSol, numTasks};
}

void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch<tempo::DefaultTime> &branch) {
    for (auto [id, val]: branch) {
        auto lit = solver.boolean.getLiteral(val, id);
        solver.set(lit);
    }

    solver.propagate();
}
