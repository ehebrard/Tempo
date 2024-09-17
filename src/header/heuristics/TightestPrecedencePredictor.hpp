/**
* @author Tim Luchterhand
* @date 16.09.24
* @brief
*/

#ifndef TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP
#define TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP

#include <vector>
#include <Iterators.hpp>

#include "util/traits.hpp"
#include "util/distance.hpp"
#include "util/IntFinity.hpp"
#include "PrecedencePredictor.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {
    /**
     * @brief Precedence predictor based on Tightest heuristic
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class TightestPrecedencePredictor : public PrecedencePredictor<TightestPrecedencePredictor<T>, T> {
    public:
        /**
         * Ctor
         * @param literals all possible search literals
         */
        explicit TightestPrecedencePredictor(std::vector<Literal<T>> literals)
                : PrecedencePredictor<TightestPrecedencePredictor<T>, T>(std::move(literals)) {}


        /**
         * Update confidences using the current state of the solver
         * @param solver Solver that represents the current state of the search
         */
        void updateConfidence(const Solver<T> &solver) {
            for(auto [lit, pos, neg] : iterators::zip(this->literals, this->massesPos, this->massesNeg)) {
                auto estPos = boundEstimation(lit, solver);
                auto estNeg = boundEstimation(~lit, solver);
                intfinity distPos = estPos.has_value() ? *estPos
                                                       : intfinity<typename decltype(estPos)::value_type>::Inf();
                intfinity distNeg = estNeg.has_value() ? *estNeg
                                                       : intfinity<typename decltype(estNeg)::value_type>::Inf();
                auto p = static_cast<float>(distPos) / static_cast<float>(distPos + distNeg);
                pos += 1 - p;
                neg += p;
            }
        }

        /**
         * Calculates prediction certainty from two evidence masses (not really theoretically motivated)
         * @param mPos positive evidence for a literal
         * @param mNeg negative evidence for a literal
         * @return confidence value in [0, 1]
         */
        static double getCertainty(double pos, double neg) {
            return std::pow(std::abs(2 * pos / (pos + neg) - 1), 0.1);
        }
    };
}

#endif //TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP
