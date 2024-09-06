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

#include "Solver.hpp"
#include "util/Options.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/serialization.hpp"
#include "Model.hpp"

/**
 * @brief represents a disjunctive resource in scheduling problems
 * @tparam T timing type
 */
template<typename T = int>
class DisjunctiveResource : public std::vector<unsigned> {
public:
    using std::vector<unsigned>::vector;

    static constexpr auto resourceCapacity() noexcept { return 1; }

    static constexpr auto getDemand(unsigned) noexcept { return 1; }
};

using SolverPtr = std::unique_ptr<tempo::Solver<>>;
using ProblemInstance = tempo::SchedulingProblemHelper<int, DisjunctiveResource<int>>;

/**
 * loads a problem instance and instantiates the solver using the given options
 * @param options options for the solver
 * @return ready to run scheduler (with default heuristics) and problem scheduling problem instance struct,
 * optionally the optimal solution and the number of tasks
 */
auto loadSchedulingProblem(
        const tempo::Options &options) -> std::tuple<SolverPtr, ProblemInstance, std::optional<int>, unsigned>;


/**
 * Puts the solver on a branch
 * @param solver solver to configure
 * @param branch branch to set the solver on
 */
void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch &branch);


#endif //TEMPO_SCHEDULING_HELPERS_HPP
