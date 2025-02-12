/**
* @author Tim Luchterhand
* @date 08.08.24
* @brief
*/

#ifndef TEMPO_VARIMPORTANCERUNNER_HPP
#define TEMPO_VARIMPORTANCERUNNER_HPP

#include <mutex>
#include <Iterators.hpp>

#include "util/Options.hpp"
#include "util/Matrix.hpp"
#include "Model.hpp"
#include "util/scheduling_helpers.hpp"
#include "util/serialization.hpp"

using Time = int;
namespace ser = tempo::serialization;

/**
 * Represents the different outcomes of a scheduler run
 */
enum class SchedulerState {
    Valid, ///< run completed successfully, the locally optimal solution was found
    Unsat, ///< problem is unsolvable
    AlreadyDecided ///< The added arc between two events is already present in or entailed by the problem definition
};

/**
 * Structure containing information of a scheduler run
 */
struct Result {
    SchedulerState state; ///< Outcome of the scheduler run
    std::optional<Time> result; ///< makespan of the solution found or null if run did not complete
};


/**
 * @brief Thread safe runner for variable importance calculation
 */
class VarImportanceRunner {
    tempo::serialization::PartialProblem problem;
    tempo::Options options;
    Time optimum;
    std::vector<bool> literalCache;
    std::vector<tempo::Literal<Time>> searchLiterals;
    std::mutex mutex;
    unsigned long totalNumDecisions = 0;
    unsigned numSearches = 0;
    bool inconsistent = false;

public:
    /**
     * Ctor
     * @param partialProblem problem instance to solve.
     * @param options scheduler options for the search
     * @param optSol optimal solution known to this problem
     */
    VarImportanceRunner(ser::PartialProblem partialProblem, tempo::Options options,
                        const ser::Solution<Time> &optSol);

    /**
     * Get the positive literals for all boolean search variables
     * @return access to the search literals
     */
    [[nodiscard]] auto getLiterals() const noexcept -> const std::vector<tempo::Literal<Time>> &;

    /**
     * Runs the solver on the problem specified ad construction while enforcing the given literal
     * @tparam Search Search function type
     * @param lit literal to enforce
     * @param searchFunction Search function object to run the actual search
     * @return result containing in run status and obtained makespan if problem + lit is SAT
     */
    template<std::invocable<tempo::Solver<Time> &, tempo::MinimizationObjective<Time>&> Search>
    auto run(tempo::Literal<Time> lit, Search &&searchFunction) -> Result {
        using enum SchedulerState;
        if (isInconsistent()) {
            throw std::runtime_error("inconsistent sub problem");
        }

        auto problemInstance = loadSchedulingProblem(options);
        auto &s = *problemInstance.solver;
        const auto &p = problemInstance.instance;
        loadBranch(s, problem.decisions);
        if (s.boolean.satisfied(lit) or s.boolean.falsified(lit)) {
            return {AlreadyDecided, {}};
        }

        {
            std::lock_guard guard(mutex);
            if (literalCache.at(lit)) {
                return {Valid, optimum};
            }
        }
        try {
            s.set(lit);
        } catch(const tempo::Failure<Time> &) {
            return {Unsat, {}};
        }

        tempo::MinimizationObjective objective(p.schedule().duration);
        s.SolutionFound.subscribe_unhandled(
            [o = optimum, durVar = p.schedule().duration, &objective](auto &solver) mutable {
                if (solver.numeric.lower(durVar) == o) {
                    objective.setDual(objective.primalBound());
                }
            });
        std::forward<Search>(searchFunction)(s, objective);
        return evaluateRun(s, p.schedule().duration);
    }

    /**
     * Get the average number of decisions over all searches
     * @return total number of decisions made divided by number of search runs
     */
    [[nodiscard]] double averageNumberOfDecisions() const noexcept;

    /**
     * Calculates the average search effort for this snapshot
     * @return average number of decisions divided by number of choice points in the problem
     */
    [[nodiscard]] double averageSearchEffort() const noexcept;

    /**
     * Whether the loaded data point is inconsistent
     * @return true if data point is inconsistent, false otherwise
     */
    [[nodiscard]] bool isInconsistent() const noexcept;

protected:


    Result evaluateRun(const tempo::Solver<Time> &solver, const tempo::NumericVar<Time> &duration);
};


#endif //TEMPO_VARIMPORTANCERUNNER_HPP
