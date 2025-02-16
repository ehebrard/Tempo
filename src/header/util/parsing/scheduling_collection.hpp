/**
* @author Tim Luchterhand
* @date 12.02.25
* @file scheduling_collection.hpp
* @brief Utils for loading arbitrary scheduling problems
*/

#ifndef SCHEDULING_COLLECTION_HPP
#define SCHEDULING_COLLECTION_HPP

#include <vector>
#include <variant>
#include <memory>
#include <optional>

#include "Model.hpp"
#include "Solver.hpp"
#include "util/Options.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/traits.hpp"
#include "util/factory_pattern.hpp"

namespace tempo {
    /**
     * @brief represents a disjunctive resource in scheduling problems
     */
    class DisjunctiveResource : public std::vector<unsigned> {
    public:
        using std::vector<unsigned>::vector;

        [[nodiscard]] constexpr auto tasks() const -> const std::vector<unsigned> & {
            return *this;
        }

        static constexpr auto resourceCapacity() noexcept { return 1; }

        static constexpr auto getDemand(unsigned) noexcept { return 1; }
    };

    /**
     * @brief represents a cumulative resource in scheduling problems
     * @tparam R resource unit type
     */
    template<concepts::scalar R = int>
    class CumulativeResource {
        std::vector<unsigned> ts;
        std::vector<R> demands;
        R capacity;

    public:
        constexpr CumulativeResource() noexcept: ts(), demands(), capacity(0) {}

        CumulativeResource(std::vector<unsigned> tasks, std::vector<R> demands,
                           R capacity) noexcept: ts(std::move(tasks)), demands(std::move(demands)),
                                                 capacity(capacity) {}

        [[nodiscard]] constexpr auto tasks() const noexcept -> const std::vector<unsigned> & {
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
    template<SchedulingResource ...Rs>
    class VariantResource : public std::variant<Rs...> {
    public:
        using std::variant<Rs...>::variant;

        [[nodiscard]] DYNAMIC_DISPATCH_VOID(tasks, const)

        [[nodiscard]] DYNAMIC_DISPATCH_VOID(resourceCapacity, const)

        [[nodiscard]] DYNAMIC_DISPATCH(getDemand, unsigned index, =, index, const)
    };

    template<resource_expression ...E>
    struct VariantResourceConstraint : public std::variant<E...> {
        using std::variant<E...>::variant;

        [[nodiscard]] DYNAMIC_DISPATCH_VOID(begin, const)

        [[nodiscard]] DYNAMIC_DISPATCH_VOID(end, const)

        [[nodiscard]] DYNAMIC_DISPATCH_VOID(getDisjunctiveLiterals, const)
    };


    using Time = int;
    using ResourceUnit = int;
    using Resource = VariantResource<DisjunctiveResource, CumulativeResource<ResourceUnit>>;
    using ResourceConstraint = VariantResourceConstraint<NoOverlapExpression<Time>, CumulativeExpression<Time>>;
    using SolverPtr = std::unique_ptr<Solver<Time>>;
    using ProblemInstance = SchedulingProblemHelper<Time, Resource>;

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
    auto loadSchedulingProblem(const Options &options) -> Problem;
}

#endif //SCHEDULING_COLLECTION_HPP
