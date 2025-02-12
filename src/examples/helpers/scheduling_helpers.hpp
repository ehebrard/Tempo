/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/

#ifndef TEMPO_SCHEDULING_HELPERS_HPP
#define TEMPO_SCHEDULING_HELPERS_HPP

#include <memory>
#include <optional>
#include <string>
#include <iostream>

#include "Model.hpp"
#include "Solver.hpp"
#include "util/Options.hpp"
#include "util/serialization.hpp"
#include "util/printing.hpp"
#include "util/Profiler.hpp"
#include "Solution.hpp"
#include "heuristics/LNS/RelaxationEvaluator.hpp"
#include "util/parsing/scheduling_collection.hpp"


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
    std::cout << "-- assumptions per run: ";
    printRange(evalPolicy.assumptionsPerRun(), std::cout) << "\n";
    std::cout << "-- runs: ";
    printRange(evalPolicy.runStatus(), std::cout) << "\n";
    std::cout << "-- fails per run: ";
    printRange(evalPolicy.failsPerRun(), std::cout) << "\n";
    std::cout << "-- search time per run: ";
    printRange(evalPolicy.searchTimesPerRun(), std::cout) << "\n";
    std::cout << "-- assumption time per run: ";
    printRange(evalPolicy.assumptionTimePerRun(), std::cout) << "\n";
    std::cout << "-- acc: ";
    printRange(evalPolicy.assumptionAccuracyPerRun(), std::cout) << "\n";
    std::cout << "-- normalized acc: ";
    printRange(evalPolicy.normalizedAssumptionAccuracyPerRun(), std::cout) << "\n";
    std::cout << "-- solution discrepancy: ";
    printRange(evalPolicy.solutionDiscrepancy(), std::cout) << "\n";
    std::cout << "-- policy run accuracy: " << evalPolicy.runAccuracy() << std::endl;
    std::cout << "-- policy assumption accuracy: " << evalPolicy.totalAssumptionAccuracy() << std::endl;
    return sw.elapsed<std::chrono::milliseconds>();
}

/**
 * @brief Helper class that maps literals to task edges
 * @details @copybrief
 */
class EdgeMapper {
    std::unique_ptr<tempo::Solver<tempo::Time>> solver;
    tempo::ProblemInstance problem;

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
    [[nodiscard]] auto getTaskEdge(tempo::Literal<tempo::Time> lit) const -> std::pair<unsigned, unsigned>;

    /**
     * @return number of tasks in the problem
     */
    [[nodiscard]] std::size_t numTasks() const noexcept;

    /**
     * Access to the solver
     * @return
     */
    auto getSolver() noexcept -> tempo::Solver<tempo::Time> &;

    /**
     * Access to the solver
     * @return
     */
    [[nodiscard]] auto getSolver() const noexcept -> const tempo::Solver<tempo::Time> &;
};


#endif //TEMPO_SCHEDULING_HELPERS_HPP
