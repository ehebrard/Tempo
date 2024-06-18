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
concept CostProvider = requires(Impl instance, var x, const S& solver) {
    { instance.getCost(x,solver) } -> std::convertible_to<double>;
};

template <typename Impl>
class RankingHeuristic : public crtp<Impl, RankingHeuristic> {
public:

    template<concepts::typed_range<var_t> Variables, BranchProvider S> requires(CostProvider<Impl, S>)
    auto bestVariable(const Variables &variables, const S &solver) const {
        auto best_var = Constant::NoVarx;
        double minCost = std::numeric_limits<double>::infinity();
        assert(not std::ranges::empty(variables));
        for (auto x: variables) {
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
    RankingHeuristic() = default;
    friend Impl;
};
}

#endif //TEMPO_RANKINGHEURISTIC_HPP
