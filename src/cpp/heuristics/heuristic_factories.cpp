/**
* @author Tim Luchterhand
* @date 09.07.24
* @brief
*/

#include <stdexcept>

#include "heuristics/heuristic_factories.hpp"

#define RETURN_NAME(NAME) return #NAME

namespace tempo::heuristics::detail {

    auto getVarHName(Options::ChoicePointHeuristics heuristic) -> std::string {
        switch (heuristic) {
            case Options::ChoicePointHeuristics::Tightest:
                RETURN_NAME(Tightest);
            case Options::ChoicePointHeuristics::WeightedDegree:
                RETURN_NAME(WeightedDegree);
            case Options::ChoicePointHeuristics::VSIDS:
                RETURN_NAME(VSIDS);
            default:
                throw std::runtime_error("unknown variable heuristic");
        }
    }

    auto getValHName(Options::PolarityHeuristic heuristic) -> std::string {
        switch (heuristic) {
            case Options::PolarityHeuristic::Tightest:
                RETURN_NAME(TightestValue);
            case Options::PolarityHeuristic::TSG:
                RETURN_NAME(TightestSolutionGuided);
            case Options::PolarityHeuristic::Random:
                RETURN_NAME(RandomBinaryValue);
            default:
                throw std::runtime_error("unknown value heuristic type");
        }
    }
}