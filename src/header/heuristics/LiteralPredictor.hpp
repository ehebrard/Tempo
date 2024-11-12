/**
* @author Tim Luchterhand
* @date 17.09.24
* @brief
*/

#ifndef TEMPO_LITERALPREDICTOR_HPP
#define TEMPO_LITERALPREDICTOR_HPP

#include <vector>
#include <concepts>

#include "util/crtp.hpp"
#include "util/traits.hpp"
#include "Literal.hpp"

namespace tempo::heuristics {

    template<typename P, typename T>
    concept LiteralWeighter = requires(const P instance, Literal<T> lit, std::size_t idx) {
        { instance.getCertainty(lit, idx) } -> concepts::scalar;
        { instance.getPolarity(lit, idx) } -> std::same_as<Literal<T>>;
    };

    /**
     * @brief Literal predictor CRTP helper class
     * @details @copybrief
     * @tparam Impl actual implementation
     * @tparam T timing type
     */
    template<typename Impl, concepts::scalar T>
    class LiteralPredictor : public crtp<Impl, LiteralPredictor, T>{
    protected:
        std::vector<Literal<T>> literals;
    public:
        /**
         * Ctor
         * @param literals all possible search literals
         */
        explicit LiteralPredictor(std::vector<Literal<T>> literals) noexcept : literals(std::move(literals)) {}

        /**
         * Gets a list of literals ordered by their prediction confidence
         * @return vector of pair(literal, confidence) ordered by confidence in descending order
         */
        auto getLiterals() const {
            static_assert(LiteralWeighter<Impl, T>);
            using namespace std::views;
            using Hax = LiteralPredictor*;
            // Dont worry! const_cast is necessary because of this issue in the standard (poison pills for ranges)
            // https://www.open-std.org/JTC1/SC22/WG21/docs/papers/2022/p2602r2.html
            auto zip = iterators::enumerate(const_cast<Hax>(this)->literals) | transform([this](auto tpl) {
                auto [idx, lit] = tpl;
                return std::make_pair(this->getImpl().getPolarity(lit, idx), this->getImpl().getCertainty(lit, idx));
            }) | common;

            std::vector ret(zip.begin(), zip.end());
            std::ranges::sort(ret, {}, [](const auto &tpl) { return -std::get<1>(tpl); });
            return ret;
        }

        /**
         * Gets a list of literals where the prediction confidence is above a given threshold
         * @param confidenceThreshold prediction certainty threshold in [0, 1]
         * @return list of literals where the prediction confidence is above confidenceThreshold
         */
        auto getLiterals(double confidenceThreshold) const -> std::vector<Literal<T>> {
            static_assert(LiteralWeighter<Impl, T>);
            std::vector<Literal<T>> ret;
            for (auto [idx, lit] : iterators::const_enumerate(literals)) {
                auto certainty = this->getImpl().getCertainty(lit, idx);
                if (certainty > confidenceThreshold) {
                    ret.emplace_back(this->getImpl().getPolarity(lit, idx));
                }
            }

            return ret;
        }

        /**
         * Number of literals
         * @return number of literals
         */
        [[nodiscard]] std::size_t numLiterals() const noexcept {
            return literals.size();
        }
    };
}

#endif //TEMPO_LITERALPREDICTOR_HPP
