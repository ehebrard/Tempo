//
// Created by tim on 16.11.22.
//

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


        auto nextVariable(Solver<T> &solver) {
            return std::visit([&solver](auto &heuristic) { return heuristic.nextVariable(solver); }, impl);
        }

    private:
        Implementations impl{};
    };
}

#endif //TEMPO_HEURISTICMANAGER_HPP

