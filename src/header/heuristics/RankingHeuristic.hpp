/**
 * @author Tim Luchterhand
 * @date 21.11.22
 */

#ifndef TEMPO_RANKINGHEURISTIC_HPP
#define TEMPO_RANKINGHEURISTIC_HPP

//#define DEBUG_HEURISTICS

#include <concepts>

#include "util/crtp.hpp"
#include "Constant.hpp"

namespace tempo::heuristics {

template<typename S>
concept BranchProvider = requires(const S solver) {
    { solver.getBranch() } -> concepts::ctyped_range<var_t>;
};

template <typename Impl, typename S>
concept CostProvider = requires(Impl instance, var_t x, const S& solver) {
    { instance.getCost(x,solver) } -> std::convertible_to<double>;
};

template <typename Impl>
class RankingHeuristic : public crtp<Impl, RankingHeuristic> {
public:

    template<concepts::typed_range<var_t> Variables, BranchProvider S> requires(CostProvider<Impl, S>)
    auto bestVariable(const Variables &variables, const S &solver) const {
        auto best_var = Constant::NoVar;
        double minCost = std::numeric_limits<double>::infinity();
        assert(not std::ranges::empty(variables));
        for (auto x: variables) {
            const auto cost = this->getImpl().getCost(x, solver);

            if (cost < minCost) {
                minCost = cost;
                best_var = x;
            }
        }

        assert(best_var != Constant::NoVar);
        return best_var;
  }

private:
    RankingHeuristic() = default;
    friend Impl;
};
}

#endif //TEMPO_RANKINGHEURISTIC_HPP
