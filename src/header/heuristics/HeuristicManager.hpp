/************************************************
 * Tempo HeuristicManager.hpp
 *
 * Copyright 2024 Tim Luchterhand  and Emmanuel Hebrard
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

#ifndef TEMPO_HEURISTICMANAGER_HPP
#define TEMPO_HEURISTICMANAGER_HPP

#include <variant>
#include <exception>
#include <utility>

#include "VSIDS.hpp"
#include "Tightest.hpp"
#include "WeightedDegree.hpp"


namespace tempo {
    template<typename T>
    class Solver;
}

/**
 * @brief namespace containing variable (choice point) selection heuristics
 */
namespace tempo::heuristics {
    /**
     * @brief Heuristic factory class that can be used to construct different heuristics and at the same time provides a
     * consistent interface to callers
     */
    template<concepts::scalar T>
    class HeuristicManager {
        using Implementations = std::variant<Tightest, VSIDS<T>, WeightedDegree<T>>;

    public:
        /**
         * Ctor: Internally constructs the heuristic inferred from the given
         * arguments
         * @param solver solver for which to create a heuristic
         * @param options options specifying the type of heuristic and further
         * config values
         * @throws std::runtime_error if an unknown heuristics type was given in
         * options
         */
        HeuristicManager(Solver<T> &solver, const Options &options) {
            //@TODO use factory pattern
            switch (options.choice_point_heuristics) {
                case Options::ChoicePointHeuristics::Tightest: {
                    impl.template emplace<Tightest>();
                    break;
                }
                case Options::ChoicePointHeuristics::VSIDS: {
                    if (options.learning) {
                        impl.template emplace<VSIDS<T>>(solver);
                    } else // closest thing if not learning
                        impl.template emplace<WeightedDegree<T>>(solver);
                    break;
                }
                case Options::ChoicePointHeuristics::WeightedDegree: {
                    impl.template emplace<WeightedDegree<T>>(solver);
                    break;
                }
                default:
                    throw std::runtime_error("unknown heuristic type");
            }
        }


        /**
         * Calls the internally stored heuristic with the given arguments
         * @param solver solver for which to select the next variable
         * @return variable choice consisting of the selected variable and its type
         */
        auto nextVariable(Solver<T> &solver) {
            return std::visit([&solver](auto &heuristic) { return heuristic.nextVariable(solver); }, impl);
        }

    private:
        Implementations impl{};
    };
}

#endif //TEMPO_HEURISTICMANAGER_HPP

