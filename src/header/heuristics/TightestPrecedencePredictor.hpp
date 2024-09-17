/**
* @author Tim Luchterhand
* @date 16.09.24
* @brief
*/

#ifndef TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP
#define TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP

#include <vector>
#include <Iterators.hpp>

#include "Literal.hpp"
#include "util/traits.hpp"
#include "util/distance.hpp"
#include "util/IntFinity.hpp"

namespace tempo::heuristics {
    template<concepts::scalar T>
    class TightestPrecedencePredictor {
        std::vector<Literal<T>> literals;
        std::vector<float> massesPos;
        std::vector<float> massesNeg;
    public:
        explicit TightestPrecedencePredictor(std::vector<Literal<T>> literals) : literals(std::move(literals)),
                                                                                 massesPos(this->literals.size(), 0),
                                                                                 massesNeg(this->literals.size(), 0) {}

        void updateConfidence(const Solver<T> &solver) {
            for(auto [lit, pos, neg] : iterators::zip(literals, massesPos, massesNeg)) {
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


        auto getLiterals(double confidenceThreshold) const -> std::vector<Literal<T>> {
            std::vector<Literal<T>> ret;
            for (auto [lit, pos, neg] : iterators::zip(literals, massesPos, massesNeg)) {
                // = |m1 - m2| / (m1 + m2)
                auto certainty = std::pow(std::abs(2 * pos / (pos + neg) - 1), 0.1);
                if (certainty > confidenceThreshold) {
                    ret.emplace_back(pos > neg ? lit : ~lit);
                }
            }

            return ret;
        }
    };
}

#endif //TEMPO_TIGHTESTPRECEDENCEPREDICTOR_HPP
