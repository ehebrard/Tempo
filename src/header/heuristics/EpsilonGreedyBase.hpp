/************************************************
 * Tempo HeuristicManager.hpp
 *
 * Copyright 2024 Tim Luchterhand 
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

#ifndef SCHEDCL_EPSILONGREEDYBASE_HPP
#define SCHEDCL_EPSILONGREEDYBASE_HPP
#include <concepts>
#include "DistanceMatrix.hpp"

namespace schedcl {
    template<typename T>
    class Scheduler;
}

namespace schedcl::heuristics {
    /**
     * @brief Requirement for a type used as parameter for EpsilonGreedyBase
     * @details @copybrief
     * Requires a member function nextChoicePoint with valid signature
     * @tparam HeuristicT
     * @tparam T
     */
    template<typename HeuristicT, typename T>
    concept Heuristic = requires(HeuristicT instance, Scheduler<T> scheduler) {
        { instance.nextChoicePoint(scheduler) } -> std::same_as<DistanceConstraint<T>>;
    };

    /**
     * @brief Base class for constructing epsilon greedy heuristics
     * @details @copybrief
     * behaves like the base heuristic with a probability of  1 - epsilon. Otherwise chooses a random choice point
     * @tparam HeuristicT type of base heuristic (e.g. VSIDS)
     */
    template<typename HeuristicT>
    class EpsilonGreedyBase {
    public:
        /**
         * CTor
         * @tparam Args argument types of base heuristic
         * @param epsilon epsilon parameter within [0, 1]
         * @param args arguments of base heurisitc
         */
        template<typename ...Args>
        constexpr explicit EpsilonGreedyBase(double epsilon, Args &&...args) :
        epsilon(static_cast<std::size_t>(epsilon * 100)), base(std::forward<Args>(args)...) {
            if (epsilon < 0 or epsilon > 1) {
                throw std::runtime_error("invalid value for epsilon");
            }
        }

        /**
         * Returns the next choice point using epsilon greedy strategy. Invokes the base heuristic with a probability
         * of 1 - epsilon. Otherwise chooses a random choice point
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to select a choice point
         * @return selected choice point
         */
        template<typename T> requires Heuristic<HeuristicT, T>
        auto nextChoicePoint(Scheduler<T> &scheduler) -> DistanceConstraint<T> {
            if (random() % 100 >= epsilon) {
                return base.nextChoicePoint(scheduler);
            }

            auto &indexSequence = scheduler.getSequence();
            const auto &choicePoints = scheduler.getChoicePoints();
            for (auto indexIt = indexSequence.rbegin(); indexIt != indexSequence.rend(); ++indexIt) {
                const auto &cp = choicePoints[*indexIt];
                if (scheduler.isGround(cp)) {
                    indexSequence.remove_back(*indexIt);
                    continue;
                }
            }

            if (indexSequence.count() == 0) {
                return DistanceConstraint<T>::none;
            }

            auto index = indexSequence.any(indexSequence.count(), random);
            indexSequence.remove_front(index);
            return choicePoints[index];
        }

    private:
        std::size_t epsilon;
        HeuristicT base;
    };
}

#endif //SCHEDCL_EPSILONGREEDYBASE_HPP
