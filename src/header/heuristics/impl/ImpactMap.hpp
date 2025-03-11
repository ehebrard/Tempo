/**
* @author Tim Luchterhand
* @date 11.03.25
* @file ImpactMap.hpp
* @brief
*/

#ifndef IMPACTMAP_HPP
#define IMPACTMAP_HPP

#include <vector>
#include <optional>

#include "util/SubscribableEvent.hpp"
#include "util/traits.hpp"
#include "Solver.hpp"

namespace tempo::heuristics::impl {
    /**
     * @brief Tracks the impact of a binary variable choice.
     * @details @copybrief
     * The impact is simply the number of additional variables fixed when setting the binary variable
     */
    class ImpactMap {
        using ImpactMeasure = float;
        std::vector<ImpactMeasure> impact;
        SubscriberHandle decisionHandle;
        SubscriberHandle conflictHandle;
        SubscriberHandle propagationCompleteHandle;
        std::optional<std::size_t> numOpenVars;
        info_t litId = Constant::NoSemantic;
    protected:
        void updateImpact(unsigned numFixed);
    public:
        ImpactMap(const ImpactMap &) = delete;
        ImpactMap(ImpactMap &&) = delete;
        ImpactMap &operator=(const ImpactMap &) = delete;
        ImpactMap &operator=(ImpactMap &&) = delete;
        ~ImpactMap() = default;

        /**
         * Ctor
         * @tparam T Timing type
         * @param solver Solver for which to track variable impact
         */
        template<concepts::scalar T>
        explicit ImpactMap(const Solver<T> &solver)
            : impact(solver.boolean.size() * 2),
              decisionHandle(solver.ChoicePoint.subscribe_handled([this](const auto &s, auto lit) {
                  numOpenVars = s.getBranch().size();
                  litId = lit.id();
              })),
              conflictHandle(solver.ConflictEncountered.subscribe_handled([this](auto &&) {
                  if (numOpenVars.has_value()) {
                      updateImpact(numOpenVars.value());
                  }

                  numOpenVars.reset();
              })),
              propagationCompleteHandle(
                  solver.PropagationCompleted.subscribe_handled([this](const auto &s) {
                      if (numOpenVars.has_value()) {
                          updateImpact(numOpenVars.value() - s.getBranch().size());
                      }

                      numOpenVars.reset();
                  })) {}


        /**
         * Get the impact of a literal choice
         * @param lit ID of the literal
         * @return impact value of lit
         */
        [[nodiscard]] auto get(info_t lit) const -> ImpactMeasure;

        /**
         * Get the impact of a literal choice
         * @tparam T timing type
         * @param lit target literal
         * @return impact value of lit
         */
        template<concepts::scalar T>
        auto get(Literal<T> lit) const -> ImpactMeasure {
            return get(lit.id());
        }

        /**
         * Get raw impact map
         * @return reference to impact map
         */
        auto getMap() const noexcept -> const std::vector<ImpactMeasure> &;
    };
}

#endif //IMPACTMAP_HPP
