/**
* @author Tim Luchterhand
* @date 24.08.24
* @brief contains util functions for measuring distances on edges
*/

#ifndef TEMPO_EDGE_DISTANCE_HPP
#define TEMPO_EDGE_DISTANCE_HPP

#include <optional>

#include "Global.hpp"
#include "Literal.hpp"
#include "Constant.hpp"
#include "DistanceConstraint.hpp"
#include "util/traits.hpp"

namespace tempo {

    /**
     * @brief Helper class that provides distance limits
     * @tparam T distance type
     */
    template<concepts::scalar T>
    struct Limits {
        /**
         * represents the infinite distance
         * @note for normal integers this is the largest possible value
         * @return infinite distance
         */
        static constexpr auto infinity() {
            if constexpr (std::numeric_limits<T>::has_infinity) {
                return std::numeric_limits<T>::infinity();
            } else {
                return std::numeric_limits<T>::max();
            }
        }

        /**
         * Represents an invalid distance
         * @note only supported for types that support NaN
         * @return invalid distance
         */
        static constexpr auto noValue() {
            static_assert(std::numeric_limits<T>::has_quiet_NaN, "Type does not support nan");
            return std::numeric_limits<T>::quiet_NaN();
        }
    };

    /**
     * Helper constant that represents invalid distance
     * @tparam T distance type
     */
    template<concepts::scalar T>
    inline constexpr auto NoDistance = Limits<T>::noValue();

    /**
     * Helper constant that represents infinite distance
     * @tparam T distance type
     */
    template<concepts::scalar T>
    inline constexpr auto InfiniteDistance = Limits<T>::infinity();

    template<typename S>
    concept distance_provider = requires(const S instance, var_t e) {
        { instance.numeric.upper(e) } -> concepts::scalar;
        { instance.numeric.lower(e) } -> concepts::scalar;
    };


    /**
     * @brief Requirement for a type that provides distances for edges between variables
     * @tparam Solver
     */
    template <typename Solver>
    concept edge_distance_provider = distance_provider<Solver> and requires(const Solver s, var_t x) {
        { s.boolean.getEdge(true, x) } -> concepts::same_template<DistanceConstraint>;
        { s.boolean.hasSemantic(x) } -> std::convertible_to<bool>;
    };

    /**
     * Computes an estimate for the length of a given edge
     * @tparam T timing type
     * @tparam S edge distance provider (usually the solver)
     * @param edge edge for which to compute length
     * @param solver solver that provides upper and lower bounds for variables
     * @return estimated length of the edge or invalid distance if edge is null-edge
     */
    template<concepts::scalar T, distance_provider S>
    auto boundEstimation(const DistanceConstraint<T> &edge, const S &solver) -> std::optional<T> {
        if (edge.isNull()) {
            return {};
        }

        return solver.numeric.upper(edge.from) - solver.numeric.lower(edge.to);
    }

    /**
     * overload of boundEstimation. Calculates an estimate of the length of an edge associated with a given literal
     * @tparam T timing type
     * @tparam S edge distance provider (usually the solver)
     * @param lit literal associated with an arc
     * @param solver solver that provides upper and lower bounds for variables
     * @return estimated length of the edge or invalid distance if edge is null-edge
     */
    template<concepts::scalar T, edge_distance_provider S>
    auto boundEstimation(Literal<T> lit, const S &solver) -> std::optional<T> {
        if (not lit.hasSemantic()) {
            return {};
        }

        return boundEstimation(solver.boolean.getEdge(lit), solver);
    }

    /**
     * overload of boundEstimation. Calculates an estimate of the length of an edge associated with a given literal
     * @tparam S edge distance provider (usually the solver)
     * @param sign sign of the literal
     * @param x corresponding variable
     * @param solver solver that provides upper and lower bounds for variables
     * @return estimated length of the edge or invalid distance if edge is null-edge
     */
    template<edge_distance_provider S>
    auto boundEstimation(bool sign, var_t x, const S &solver) {
        using T = decltype(boundEstimation(solver.boolean.getEdge(sign, x), solver));
        if (not solver.boolean.hasSemantic(x)) {
            T empty;
            return empty;
        }

        auto edge = solver.boolean.getEdge(sign, x);
        return boundEstimation(edge, solver);
    }
}

#endif //TEMPO_EDGE_DISTANCE_HPP
