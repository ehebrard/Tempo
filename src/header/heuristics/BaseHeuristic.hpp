/**
 * @author Tim Luchterhand
 * @date 21.11.22
 */

#ifndef TEMPO_BASEHEURISTIC_HPP
#define TEMPO_BASEHEURISTIC_HPP

//#define DEBUG_HEURISTICS

#include <concepts>

#include "util/crtp.hpp"
#include "Constant.hpp"

namespace tempo::heuristics {

template<typename S>
concept BranchProvider = requires(const S solver) {
    { solver.getBranch() } -> concepts::ctyped_range<var_t>;
};

/**
 * @brief Requirement for a class that derives from BaseHeuristic
 * @details @copybrief
 * Requires a member function named getCost with valid signature
 * @tparam Impl
 */
template <typename Impl, typename S>
concept HeuristicImplementation = requires(Impl instance, var x, const S& solver) {
    { instance.getCost(x,solver) } -> std::convertible_to<double>;
};

    /**
     * @brief CRTP base class used to build heuristics that select a choice point by minimizing a cost function
     * @tparam Impl type of concrete Implementation
     */
template <typename Impl>
class BaseHeuristic : public crtp<Impl, BaseHeuristic> {
public:
  /**
   * @brief Static polymorphic interface for the caller of the heuristic.
   * Selects a choice point from a scheduler by minimizing a cost function.
   * @details @copybrief
   * Also removes ground instances from the index sequence of the scheduler
   * @tparam T type of scheduler
   * @param solver solver for which to generate a choice point
   * @return selected choice point or DistanceConstraint::none if no further
   * choices can be made
   */
    template<BranchProvider S> requires(HeuristicImplementation<Impl, S>)
    auto nextChoicePoint(const S &solver) {
        auto best_var = Constant::NoVarx;
        const auto &indexSequence = solver.getBranch();
        double minCost = std::numeric_limits<double>::infinity();
        assert(not indexSequence.empty());
        for (auto x: indexSequence) {
            const auto cost = this->getImpl().getCost(x, solver);

#ifdef DEBUG_HEURISTICS_CHOICE
            std::cout << scheduler.getEdge(POS(x)) << "<>" << scheduler.getEdge(NEG(x)) << ": " << cost;
#endif

            if (cost < minCost) {
                minCost = cost;
                best_var = x;

#ifdef DEBUG_HEURISTICS_CHOICE
                std::cout << "*";
#endif
            }

#ifdef DEBUG_HEURISTICS_CHOICE
            std::cout << std::endl;
#endif
        }

        assert(best_var != Constant::NoVarx);
        return best_var;
  }

private:
    BaseHeuristic() = default;
    friend Impl;
};
}

#endif //TEMPO_BASEHEURISTIC_HPP
