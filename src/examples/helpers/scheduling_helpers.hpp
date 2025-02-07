/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/

#ifndef TEMPO_SCHEDULING_HELPERS_HPP
#define TEMPO_SCHEDULING_HELPERS_HPP

#include <vector>
#include <memory>
#include <tuple>
#include <optional>
#include <ranges>
#include <variant>

#include "Solver.hpp"
#include "util/Options.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/serialization.hpp"
#include "util/factory_pattern.hpp"
#include "util/printing.hpp"
#include "Solution.hpp"
#include "Model.hpp"
#include "heuristics/LNS/RelaxationEvaluator.hpp"

/**
 * @brief represents a disjunctive resource in scheduling problems
 */
class DisjunctiveResource : public std::vector<unsigned> {
public:
    using std::vector<unsigned>::vector;

    [[nodiscard]] constexpr auto tasks() const -> const std::vector<unsigned>& {
        return *this;
    }

    static constexpr auto resourceCapacity() noexcept { return 1; }

    static constexpr auto getDemand(unsigned) noexcept { return 1; }
};

/**
 * @brief represents a cumulative resource in scheduling problems
 * @tparam R resource unit type
 */
template<tempo::concepts::scalar R = int>
class CumulativeResource {
    std::vector<unsigned> ts;
    std::vector<R> demands;
    R capacity;
public:
    constexpr CumulativeResource() noexcept: ts(), demands(), capacity(0) {}

    CumulativeResource(std::vector<unsigned> tasks, std::vector<R> demands, R capacity) noexcept:
            ts(std::move(tasks)), demands(std::move(demands)), capacity(capacity) {}

    [[nodiscard]] constexpr auto tasks() const noexcept -> const std::vector<unsigned> &{
        return ts;
    }

    [[nodiscard]] constexpr auto resourceCapacity() const noexcept {
        return capacity;
    }

    [[nodiscard]] constexpr auto getDemand(std::size_t taskIndex) const noexcept {
        return demands[taskIndex];
    }
};

/**
 * @brief Variant wrapper around different types of resources
 * @tparam Rs resource types
 */
template<tempo::SchedulingResource ...Rs>
class VariantResource : public std::variant<Rs...> {
public:
    using std::variant<Rs...>::variant;

    [[nodiscard]] DYNAMIC_DISPATCH_VOID(tasks, const)

    [[nodiscard]] DYNAMIC_DISPATCH_VOID(resourceCapacity, const)

    [[nodiscard]] DYNAMIC_DISPATCH(getDemand, unsigned index, =, index, const)
};

template<tempo::resource_expression ...E>
struct VariantResourceConstraint : public std::variant<E...> {
    using std::variant<E...>::variant;

    [[nodiscard]] DYNAMIC_DISPATCH_VOID(begin, const)

    [[nodiscard]] DYNAMIC_DISPATCH_VOID(end, const)

    [[nodiscard]] DYNAMIC_DISPATCH_VOID(getDisjunctiveLiterals, const)
};


using Time = int;
using ResourceUnit = int;
using Resource = VariantResource<DisjunctiveResource, CumulativeResource<ResourceUnit>>;
using ResourceConstraint = VariantResourceConstraint<tempo::NoOverlapExpression<Time>,
                                                     tempo::CumulativeExpression<Time>>;
using SolverPtr = std::unique_ptr<tempo::Solver<Time>>;
using ProblemInstance = tempo::SchedulingProblemHelper<Time, Resource>;

struct Problem {
    SolverPtr solver;
    ProblemInstance instance;
    std::vector<ResourceConstraint> constraints;
    std::optional<Time> optimalSolution;
    Time upperBound;
    unsigned numTasks;
};

/**
 * loads a problem instance and instantiates the solver using the given options
 * @param options options for the solver
 * @return ready to run scheduler (with default heuristics) and problem scheduling problem instance struct,
 * optionally the optimal solution and the number of tasks
 */
auto loadSchedulingProblem(const tempo::Options &options) -> Problem;


/**
 * Puts the solver on a branch
 * @param solver solver to configure
 * @param branch branch to set the solver on
 */
void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch &branch);

/**
 * Converts a serialization solution
 * @tparam T timing type
 * @param solution solution
 * @param options solver options
 * @return converted solution
 */
template<tempo::concepts::scalar T>
auto toSolution(const tempo::serialization::Solution<T> &solution,
                const tempo::Options &options) -> tempo::Solution<T> {
    auto problem = loadSchedulingProblem(options);
    loadBranch(*problem.solver, solution.decisions);
    problem.solver->saveSolution();
    return tempo::Solution(*problem.solver);
}

/**
 * @tparam Policy LNS relaxation policy type
 * @tparam T timing type
 * @param policy LNS relaxation policy
 * @param solver solver instance
 * @param objective objective variable
 * @param stats whether to record policy stats
 * @param solPath path to optimal solution (ignored if empty)
 * @return elapsed time in ms
 */
template<tempo::lns::relaxation_policy Policy, tempo::concepts::scalar T>
auto runLNS(Policy &&policy, tempo::Solver<T> &solver,
            tempo::MinimizationObjective<T> objective, bool stats, const std::string &solPath) {
    namespace ser = tempo::serialization;
    using namespace tempo;
    util::StopWatch sw;
    if (not stats) {
        solver.largeNeighborhoodSearch(objective, std::forward<Policy>(policy));
        return sw.elapsed<std::chrono::milliseconds>();
    }

    std::optional<Solution<Time>> solution;
    if (not solPath.empty()) {
        const auto sol = ser::deserializeFromFile<ser::Solution<int>>(solPath);
        solution = toSolution(sol, solver.getOptions());
    }

    auto evalPolicy = lns::make_evaluator(std::forward<Policy>(policy), std::move(solution));
    sw.start();
    solver.largeNeighborhoodSearch(objective, evalPolicy);
    std::cout << "-- policy run accuracy: " << evalPolicy.runAccuracy() << std::endl;
    std::cout << "-- policy assumption accuracy: " << evalPolicy.totalAssumptionAccuracy() << std::endl;
    std::cout << "-- runs: ";
    printRange(evalPolicy.runStatus(), std::cout) << "\n";
    std::cout << "-- acc: ";
    printRange(evalPolicy.assumptionAccuracyPerRun(), std::cout) << "\n";
    std::cout << "-- normalized acc: ";
    printRange(evalPolicy.normalizedAssumptionAccuracyPerRun(), std::cout) << "\n";
    std::cout << "-- assumptions per run: ";
    printRange(evalPolicy.assumptionsPerRun(), std::cout) << "\n";
    std::cout << "-- solution discrepancy: ";
    printRange(evalPolicy.solutionDiscrepancy(), std::cout) << "\n";
    std::cout << "-- fails per run: ";
    printRange(evalPolicy.failsPerRun(), std::cout) << "\n";
    std::cout << "-- search time per run: ";
    printRange(evalPolicy.searchTimesPerRun(), std::cout) << "\n";
    return sw.elapsed<std::chrono::milliseconds>();
}

/**
 * @brief Helper class that maps literals to task edges
 * @details @copybrief
 */
class EdgeMapper {
    std::unique_ptr<tempo::Solver<Time>> solver;
    ProblemInstance problem;

public:
    /**
     * Ctor
     * @param options options for solver
     */
    explicit EdgeMapper(const tempo::Options &options);

    /**
     * @brief
     * @param lit literal to map
     * @return corresponding TASK edge
     */
    [[nodiscard]] auto getTaskEdge(tempo::Literal<Time> lit) const -> std::pair<unsigned, unsigned>;

    /**
     * @return number of tasks in the problem
     */
    [[nodiscard]] std::size_t numTasks() const noexcept;

    /**
     * Access to the solver
     * @return
     */
    auto getSolver() noexcept -> tempo::Solver<Time> &;

    /**
     * Access to the solver
     * @return
     */
    [[nodiscard]] auto getSolver() const noexcept -> const tempo::Solver<Time> &;
};


#endif //TEMPO_SCHEDULING_HELPERS_HPP
