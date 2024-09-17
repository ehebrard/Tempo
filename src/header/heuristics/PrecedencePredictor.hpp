/**
* @author Tim Luchterhand
* @date 17.09.24
* @brief
*/

#ifndef TEMPO_PRECEDENCEPREDICTOR_HPP
#define TEMPO_PRECEDENCEPREDICTOR_HPP

#include <vector>
#include <concepts>

#include "util/crtp.hpp"
#include "util/traits.hpp"
#include "Literal.hpp"

namespace tempo::heuristics {

    template<typename Impl>
    concept precedence_predictor_impl = requires(const Impl impl, double d) {
        { impl.getCertainty(d, d)} -> std::convertible_to<double>;
    };

    template<typename Impl, concepts::scalar T>
    class PrecedencePredictor : public crtp<Impl, PrecedencePredictor, T>{
    protected:
        std::vector<Literal<T>> literals;
        std::vector<double> massesPos;
        std::vector<double> massesNeg;
    public:
        explicit PrecedencePredictor(std::vector<Literal<T>> literals) : literals(std::move(literals)),
                                                                         massesPos(this->literals.size(), 0),
                                                                         massesNeg(this->literals.size(), 0) {}

        auto getLiterals() const -> std::vector<std::pair<Literal<T>, double>> {
            using namespace std::views;
            using Hax = PrecedencePredictor*;
            // Dont worry! const_cast is necessary because of this issue in the standard (poison pills for ranges)
            // https://www.open-std.org/JTC1/SC22/WG21/docs/papers/2022/p2602r2.html
            auto zip = iterators::zip(const_cast<Hax>(this)->literals, const_cast<Hax>(this)->massesPos,
                                      const_cast<Hax>(this)->massesNeg) | transform([this](auto tpl) {
                auto [lit, pos, neg] = tpl;
                return std::make_pair(pos > neg ? lit : ~lit, this->getImpl().getCertainty(pos, neg));
            });

            std::vector ret(zip.begin(), zip.end());
            std::ranges::sort(ret, {}, [](const auto &tpl) { return -std::get<1>(tpl); });
            return ret;
        }

        auto getLiterals(double confidenceThreshold) const -> std::vector<Literal<T>> {
            std::vector<Literal<T>> ret;
            for (auto [lit, pos, neg] : iterators::zip(literals, massesPos, massesNeg)) {
                // = |m1 - m2| / (m1 + m2)
                auto certainty = this->getImpl().getCertainty(pos, neg);
                if (certainty > confidenceThreshold) {
                    ret.emplace_back(pos > neg ? lit : ~lit);
                }
            }

            return ret;
        }
    };
}

#endif //TEMPO_PRECEDENCEPREDICTOR_HPP
