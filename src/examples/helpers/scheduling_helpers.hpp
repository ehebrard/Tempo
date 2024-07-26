/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief
*/

#ifndef TEMPO_SCHEDULING_HELPERS_HPP
#define TEMPO_SCHEDULING_HELPERS_HPP

#include <vector>
#include <memory>
#include <utility>

#include "Solver.hpp"
#include "util/Options.hpp"
#include "util/SchedulingProblemHelper.hpp"
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
 * @return ready to run scheduler (with default heuristics) and problem scheduling problem instance struct
 */
auto loadSchedulingProblem(const tempo::Options &options) -> std::pair<SolverPtr, ProblemInstance>;

#endif //TEMPO_SCHEDULING_HELPERS_HPP
