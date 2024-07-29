/**
* @author Tim Luchterhand
* @date 24.07.24
* @brief
*/

#ifndef TEMPO_TRACEWATCHER_HPP
#define TEMPO_TRACEWATCHER_HPP

#include <vector>
#include <Iterators.hpp>
#include <cassert>
#include <ranges>
#include <algorithm>

#include "util/traits.hpp"
#include "Model.hpp"

namespace tempo {

    /**
     * @brief Class that memorizes choices taken on the path to a solution.
     * @details @copybrief Use it to determine if a partial solution is a subset of a previously found solution
     * @note only considers binary literals
     */
    class TraceWatcher {
        std::vector<bool> varPolarity;
        bool onTrack;
        var_t offset;
    public:

        /**
         * CTor
         * @param vars range containing ids of all search variables
         */
        template<concepts::ctyped_range<var_t> Vars>
        requires(std::ranges::sized_range<Vars>)
        explicit TraceWatcher(const Vars &vars): varPolarity(std::ranges::size(vars), false), onTrack(false),
                                                 offset(std::ranges::min(vars)) {
            if (std::ranges::max(vars) - offset + 1 != varPolarity.size()) {
                throw std::runtime_error("expected continuous range of search variables");
            }
        }

        /**
         * registers a complete solution
         * @tparam TF truth function type
         * @param truthFunction function object that maps variable ids to boolean values
         */
        template<concepts::callable_r<bool, var_t> TF>
        void registerSolution(TF &&truthFunction) {
            for (auto [var, val] : iterators::enumerate(varPolarity, offset)) {
                val = std::forward<TF>(truthFunction)(var);
            }

            onTrack = true;
        }

        /**
         * Checks if still on track to last solution and sets the internal on track flag to false if this is not the
         * case. Also returns a list with variables in conflict with the previous solution along with their
         * corresponding truth values
         * @tparam AF alignment function
         * @param isAligned function object that checks if the truth value of a variable aligns with a given boolean.
         * For example: undefined aligns with true and false but false only aligns with false.
         * @return Vector of pairs of conflicting variables along with their corresponding truth value from the last
         * solution
         */
        template<concepts::callable_r<bool, var_t, bool> AF>
        auto updateOnTrack(AF &&isAligned) -> std::vector<std::pair<var_t, bool>> {
            std::vector<std::pair<var_t, bool>> ret;
            for (auto [var, val] : iterators::const_enumerate(varPolarity, offset)) {
                if (not std::forward<AF>(isAligned)(var, val)) {
                    ret.emplace_back(var, val);
                }
            }

            onTrack = ret.empty();
            return ret;
        }

        /**
         * Registers a decision
         * @tparam T timing type
         * @param decision decision literal
         */
        template<concepts::scalar T>
        void step(Literal<T> decision) noexcept {
            assert(decision.isBoolean());
            onTrack &= varPolarity[decision.variable() - offset] == decision.sign();
        }

        /**
         * Manually sets the on track flag
         * @param truthVal target value
         */
        void setOnTrack(bool truthVal) noexcept;

        /**
         * Whether still on track to last solution
         * @return true if still on track to last solution, false otherwise
         */
        [[nodiscard]] bool isOnTrack() const noexcept;

        /**
         * Get the last registered solution
         * @return Vector with truth value for each variable (indexed by variable id)
         */
        [[nodiscard]] auto getLastSolution() const noexcept -> const std::vector<bool>&;

        [[nodiscard]] var_t getOffset() const noexcept;

    };
}

#endif //TEMPO_TRACEWATCHER_HPP
