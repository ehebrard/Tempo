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
#include "Model.hpp"

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
    [[nodiscard]] constexpr decltype(auto) tasks() const {
        return std::visit([](const auto &r) -> decltype(auto) { return r.tasks(); }, *this);
    }

    [[nodiscard]] constexpr auto resourceCapacity() const {
        return std::visit([](const auto &r) { return r.resourceCapacity(); }, *this);
    }

    [[nodiscard]] constexpr auto getDemand(unsigned taskIndex) const {
        return std::visit([taskIndex](const auto &r) { return r.getDemand(taskIndex); }, *this);
    }
};

template<tempo::resource_expression ...E>
struct VariantResourceConstraint : public std::variant<E...> {
    using std::variant<E...>::variant;
    constexpr auto begin() const noexcept {
        return std::visit([](const auto &e) { return e.begin(); }, *this);
    }

    constexpr auto end() const noexcept {
        return std::visit([](const auto &e) { return e.end(); }, *this);
    }

    constexpr auto begDisjunct() const noexcept {
        return std::visit([](const auto &e) { return e.begDisjunct(); }, *this);
    }

    constexpr auto endDisjunct() const noexcept {
        return std::visit([](const auto &e) { return e.endDisjunct(); }, *this);
    }
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


#endif //TEMPO_SCHEDULING_HELPERS_HPP
