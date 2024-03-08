//
// Created by tluchterha on 24/11/22.
//

#ifndef SCHEDCL_DECAYINGACTIVITYMAP_HPP
#define SCHEDCL_DECAYINGACTIVITYMAP_HPP
#include "heuristics/impl/ActivityMap.hpp"
#include "DistanceMatrix.hpp"

namespace schedcl::heuristics::impl {
    /**
     * @brief Records activity for literals. Activity of inactive literals decays which each update
     */
    class DecayingActivityMap : public ActivityMap {
    public:
        /**
         * @copydoc ActivityMap::ActivityMap
         * @param decay decay value. Each entry is multiplied by this value, each time the activity is updated
         */
        template<typename T>
        explicit DecayingActivityMap(const Scheduler<T> &scheduler, double decay) :
        ActivityMap(scheduler), decay(decay) {
            if (decay <= 0 or decay > 1) {
                throw std::runtime_error("invalid decay value");
            }
        }

        /**
         * Updates the activity map of the variables involved in the given clause, also applies the decay to each value
         * @tparam Clause
         * @param clause
         */
        template<typename Clause>
        void update(const Clause &clause) noexcept {

          //            std::cout << "update activity\n";

          static_assert(traits::is_iterable_v<Clause>,
                        "Expected iterable for Clause type");
          static_assert(
              traits::is_same_template_v<DistanceConstraint,
                                         traits::value_type_t<Clause>>,
              "clause must contain DistanceConstraints");
          bool normalize = false;
          for (const auto &arc : clause) {
            // ignore clauses that are not within the choice points node
            // interval
            if (this->contains(arc)) {
              auto &activity = this->get(arc);
              activity += increment;
              if (activity > maxActivity) {
                normalize = true;
              }
            }
          }

            // protect against overflow
            if (normalize) {
                this->for_each([i = increment](auto &val) { val /= i; });
                increment = 1;
            }

            increment /= decay;

            //            this->display(std::cout);
            //            std::cout << std::endl;
        }

    private:
        static constexpr double maxActivity = 1e10;
        double decay;
        double increment = 1;
    };
}

#endif //SCHEDCL_DECAYINGACTIVITYMAP_HPP
