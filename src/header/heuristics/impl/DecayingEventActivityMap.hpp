//
// Created by tluchterha on 24/11/22.
//

#ifndef TEMPO_DECAYINGEVENTACTIVITYMAP_HPP
#define TEMPO_DECAYINGEVENTACTIVITYMAP_HPP
#include "heuristics/impl/EventActivityMap.hpp"
//#include "DistanceMatrix.hpp"

namespace tempo::heuristics::impl {
    /**
     * @brief Records activity for literals. Activity of inactive literals decays which each update
     */
template<typename T>
    class DecayingEventActivityMap : public EventActivityMap<T> {
    public:
        /**
         * @copydoc ActivityMap::ActivityMap
         * @param decay decay value. Each entry is multiplied by this value, each time the activity is updated
         */
//        template<typename T>
        explicit DecayingEventActivityMap(const Scheduler<T> &scheduler, double decay) :
        EventActivityMap<T>(scheduler), decay(decay) {
            if (decay <= 0 or decay > 1) {
                throw std::runtime_error("invalid decay value");
            }
        }
        
        
        bool incrementActivity(const event x) noexcept {
            
#ifdef DEBUG_HEURISTICS
                std::cout << "  ++a[" << prettyEvent(x) << "]" << std::endl;
#endif
            
            EventActivityMap<T>::activity[x] += increment;
            return (EventActivityMap<T>::activity[x] > maxActivity);
        }
            
        
        /**
         * Updates the activity map of the variables involved in the given clause, also applies the decay to each value
         * @tparam Clause
         * @param clause
         */
        template<class Iterable>
        void update(Iterable &clause) noexcept {
//            static_assert(traits::is_iterable_v<Clause>, "Expected iterable for Clause type");
//            static_assert(traits::is_same_template_v<DistanceConstraint, traits::value_type_t<Clause>>,
//                          "clause must contain DistanceConstraints");
            bool normalize = false;
            
//            std::cout << "update (" << clause.size() << ")\n";
            
//            conflict.explain(NoLit, clause);
            for (const auto gl : clause) {
                
#ifdef DEBUG_HEURISTICS
                std::cout << " increment activity because of " << EventActivityMap<T>::sched.prettyLiteral(gl) << std::endl;
#endif
                
                if(LTYPE(gl) == BOUND_LIT) {
                    normalize |= incrementActivity(EVENT(EventActivityMap<T>::sched.getBoundLiteral(FROM_GEN(gl))));
                } else {
                    DistanceConstraint<T> bc{EventActivityMap<T>::sched.getEdge(FROM_GEN(gl))};
                    normalize |= incrementActivity(bc.from);
                    normalize |= incrementActivity(bc.to);
                }
                
            }
            
#ifdef DEBUG_HEURISTICS
            for(event x{0}; x<static_cast<event>(EventActivityMap<T>::activity.size()); ++x) {
                std::cout << std::setw(6) << prettyEvent(x);
            }
            std::cout << std::endl;
            for(event x{0}; x<static_cast<event>(EventActivityMap<T>::activity.size()); ++x) {
                std::cout << std::setw(6) << EventActivityMap<T>::activity[x];
            }
            std::cout << std::endl;
 #endif

            // protect against overflow
            if (normalize) {
                
#ifdef DEBUG_HEURISTICS
                std::cout << "\nnormalize (" << increment << ")\n" ;
#endif

                auto u{
                    *(std::max_element(EventActivityMap<T>::activity.begin(),
                                       EventActivityMap<T>::activity.end()))};
                auto l{
                    *(std::min_element(EventActivityMap<T>::activity.begin(),
                                       EventActivityMap<T>::activity.end()))};

                this->for_each([lb = l, gap = u - l](auto &val) {
                  val = (val - lb) / gap * baseGap + baseIncrement;
                });
                increment = baseIncrement;
            }

            increment /= decay;
            
//            clause.clear();
//            this->display(std::cout);
//            std::cout << std::endl;
        }

    private:
//        std::vector<lit> clause;
        static constexpr double maxActivity = 1e12;
        static constexpr double baseIncrement = 1;
        static constexpr double baseGap = 999;
        double decay;
        double increment = baseIncrement;
    };
}

#endif //SCHEDCL_DECAYINGEVENTACTIVITYMAP_HPP

// x - y <= a
// x <= b
//
