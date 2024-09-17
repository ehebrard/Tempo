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
    template<concepts::scalar T>
    class TightestPrecedencePredictor : public PrecedencePredictor<TightestPrecedencePredictor<T>, T> {
    public:
        explicit TightestPrecedencePredictor(std::vector<Literal<T>> literals)
                : PrecedencePredictor<TightestPrecedencePredictor<T>, T>(std::move(literals)) {}

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

        static double getCertainty(double pos, double neg) {
            return std::pow(std::abs(2 * pos / (pos + neg) - 1), 0.1);
        }
    };
}

#endif //TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP
