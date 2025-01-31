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


static auto parseDisjunctive(const tempo::Options &options, tempo::Solver<Time> &solver, tempo::Interval<Time> schedule)
        -> std::optional<std::pair<ProblemInstance , std::optional<Time>>> {
    std::optional<Time> optSol;
    std::vector<tempo::Interval<Time>> tasks;
    std::vector<DisjunctiveResource> resources;
    std::vector<tempo::DistanceConstraint<Time>> precedences;
    if (options.input_format == "osp") {
        optSol = osp::parse(options.instance_file, solver, schedule, tasks, resources);
    } else if (options.input_format == "jsp") {
        jsp::parse(options.instance_file, solver, schedule, tasks, resources, &precedences);
    } else if (options.input_format == "jstl") {
        jstl::parse(options.instance_file, solver, schedule, tasks, resources, &precedences);
    } else {
        return {};
    }

    auto tview = resources | std::views::transform([](auto &r) {return Resource(std::move(r));});
    std::vector<Resource> res(tview.begin(), tview.end());
    return {{ProblemInstance(std::move(tasks), std::move(res), std::move(precedences), schedule), optSol}};
}

static auto parseCumulative(const tempo::Options &options, tempo::Solver<Time> &solver, tempo::Interval<Time> schedule)
        -> std::optional<std::pair<ProblemInstance, std::optional<Time>>> {
    std::vector<tempo::Interval<Time>> tasks;
    std::vector<std::vector<std::size_t>> resources; // elems are the resource requirements of a task
    std::vector<std::vector<int>> demands; // elems are the resource demands of a task
    std::vector<std::pair<int, int>> p;
    std::vector<std::vector<int>> g;
    std::vector<ResourceUnit> capacities; // resource capacities
    std::vector<tempo::DistanceConstraint<Time>> precedences;
    if (options.input_format == "rcpsp") {
        rcpsp::parse(options.instance_file, solver, schedule, tasks, resources, demands, capacities, p, g, &precedences);
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

static auto postDisjunctive(tempo::Solver<Time> &solver,
                            const ProblemInstance &problem) -> std::vector<ResourceConstraint> {
    const auto &tasks = problem.tasks();
    std::vector<ResourceConstraint> constraints;
    constraints.reserve(problem.resources().size());
    for (const auto &consumingTasks: problem.resources()) {
        auto constraint = tempo::NoOverlap(problem.schedule(),
                                           consumingTasks.tasks() |
                                           std::views::transform([&tasks](auto taskId) { return tasks.at(taskId); }),
                                           std::vector<std::vector<Time>>{});
        solver.post(constraint);
        constraints.emplace_back(constraint);
    }

    return constraints;
}

static auto postCumulative(tempo::Solver<Time> &solver, const ProblemInstance &problem) -> std::vector<ResourceConstraint> {
    using namespace std::views;
    using namespace tempo;
    std::vector<ResourceConstraint> constraints;
    constraints.reserve(problem.resources().size());
    for (const auto &resource: problem.resources()) {
        auto capacity = solver.newConstant(resource.resourceCapacity());
        auto tView = resource.tasks() | transform([&problem](auto taskId) { return problem.tasks().at(taskId); });
        auto demView = iota(0u, static_cast<unsigned>(resource.tasks().size())) |
                       transform([&resource, &solver](auto idx) {
                           return solver.newConstant(resource.getDemand(idx));
                       });
        std::vector<Interval<Time>> tasks;
        tasks.reserve(tView.size());
        std::vector<NumericVar<ResourceUnit>> demands;
        demands.reserve(demView.size());
        for (auto [t, d]: iterators::zip(tView, demView)) {
            tasks.emplace_back(t);
            demands.emplace_back(d);
        }

        auto constraint = Cumulative(problem.schedule(), capacity, tasks, demands);
        solver.post(constraint);
        constraints.emplace_back(constraint);
    }

    return constraints;
}

auto loadSchedulingProblem(const tempo::Options &options)
-> Problem {
    using namespace tempo;
    using namespace std::views;
    auto solver = std::make_unique<Solver<>>(options);
    const auto schedule{solver->newInterval(0, Constant::Infinity<Time>, 0, 0, 0, Constant::Infinity<Time>)};
    auto res = parseDisjunctive(options, *solver, schedule);
    std::vector<ResourceConstraint> constraints;
    if (res.has_value()) {
        auto &[problem, optSol] = *res;
        constraints = postDisjunctive(*solver, problem);
    } else if ((res = parseCumulative(options, *solver, schedule)).has_value()) {
        auto [problem, _] = *res;
        constraints = postCumulative(*solver, problem);
    } else {
        std::cerr << "problem type " << options.input_format << " is not (yet) supported" << std::endl;
        std::exit(1);
    }

    auto &[problem, optSol] = *res;
    Time ub = options.ub;
    if (ub == -1) {
        ub = 0;
        for (const auto &t : problem.tasks()) {
            auto maxDur = t.maxDuration(*solver);
            if (maxDur == Constant::Infinity<Time>) {
                ub = maxDur;
                break;
            }

            ub += maxDur;
        }
    }

    solver->set(schedule.end.before(ub));
    auto numTasks = static_cast<unsigned>(problem.tasks().size());
    return {.solver = std::move(solver), .instance = std::move(problem), .constraints = std::move(
            constraints), .optimalSolution = optSol, .upperBound = ub, .numTasks = numTasks};
}

void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch &branch) {
    for (auto [id, val] : branch) {
        auto lit = solver.boolean.getLiteral(val, id);
        solver.set(lit);
    }

    solver.propagate();
}


EdgeMapper::EdgeMapper(const tempo::Options &options) {
    auto res = loadSchedulingProblem(options);
    solver = std::move(res.solver);
    problem = std::move(res.instance);
}

auto EdgeMapper::getTaskEdge(tempo::Literal<Time> lit) const -> std::pair<unsigned, unsigned> {
    if (not lit.hasSemantic()) {
        throw std::runtime_error("cannot get edge without semantic");
    }

    auto edge = solver->boolean.getEdge(lit);
    const auto &mapping = problem.getMapping();
    if (not mapping.contains(edge.from) or not mapping.contains(edge.to)) {
        throw std::runtime_error("edge does not correspond to task - task edge");
    }

    return {mapping(edge.from), mapping(edge.to)};
}

std::size_t EdgeMapper::numTasks() const noexcept {
    return problem.tasks().size();
}

auto EdgeMapper::getSolver() noexcept -> tempo::Solver<Time> & {
    return *solver;
}

auto EdgeMapper::getSolver() const noexcept -> const tempo::Solver<Time> & {
    return *solver;
}
