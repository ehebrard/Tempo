/**
* @author Tim Luchterhand
* @date 08.08.24
* @brief
*/

#ifndef TEMPO_VARIMPORTANCERUNNER_HPP
#define TEMPO_VARIMPORTANCERUNNER_HPP

#include "util/Options.hpp"
#include "util/Matrix.hpp"
#include "Model.hpp"
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


class VarImportanceRunner {
    tempo::serialization::PartialProblem problem;
    tempo::Options options;
    Time optimum;
    std::vector<bool> literalCache;
    std::vector<tempo::Literal<Time>> searchLiterals;
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
     * @param lit literal to enforce
     * @return result containing in run status and obtained makespan if problem + lit is SAT
     */
    auto run(tempo::Literal<Time> lit) -> Result;

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

};


#endif //TEMPO_VARIMPORTANCERUNNER_HPP
