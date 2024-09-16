/************************************************
 * Tempo RankingHeuristic.hpp
 *
 * Copyright 2024 Tim Luchterhand and Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#ifndef TEMPO_RANKINGHEURISTIC_HPP
#define TEMPO_RANKINGHEURISTIC_HPP


#include <concepts>
#include <ranges>

#include "heuristic_interface.hpp"
#include "util/crtp.hpp"
#include "Constant.hpp"

namespace tempo::heuristics {

template<typename S>
concept BranchProvider = requires(const S solver) {
    { solver.getBranch() } -> concepts::ctyped_range<var_t>;
};

template <typename Impl, typename S>
concept PartialOrder = requires(Impl instance, var_t x, const S& solver) {
    { instance.chooseBest(x, x, solver) } -> std::same_as<var_t>;
};

template <typename Impl>
class RankingHeuristic : public crtp<Impl, RankingHeuristic> {
public:

    template<concepts::typed_range<var_t> Variables, typename S> requires(PartialOrder<Impl, S>)
    auto bestVariable(const Variables &variables, const S &solver) const {
        assert(not std::ranges::empty(variables));
        auto best_var = *std::begin(variables);
        for (auto x: variables | std::views::drop(1)) {
            best_var = this->getImpl().chooseBest(best_var, x, solver);
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
