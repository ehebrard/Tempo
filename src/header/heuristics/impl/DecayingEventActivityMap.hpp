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
        explicit DecayingEventActivityMap(Scheduler<T> &scheduler, double decay) :
        EventActivityMap<T>(scheduler), decay(decay) {
            if (decay <= 0 or decay > 1) {
                throw std::runtime_error("invalid decay value");
            }
        }
        
        explicit DecayingEventActivityMap(Solver<T> &solver, double decay) :
        EventActivityMap<T>(solver), decay(decay) {
            if (decay <= 0 or decay > 1) {
                throw std::runtime_error("invalid decay value");
            }
        }
        
        
        bool incrementActivity(const event x) noexcept {
            
#ifdef DEBUG_HEURISTICS
                std::cout << "  ++a[" << prettyEvent(x) << "]" << std::endl;
#endif
            
            EventActivityMap<T>::numeric_activity[x] += increment;
            return (EventActivityMap<T>::numeric_activity[x] > maxActivity);
        }
        
        bool incrementActivity(const Literal<T> l, const Solver<T>& solver) noexcept {
            
//            std::cout << "increment activity of " << l << std::endl;
            
            bool need_rescaling{false};
            if(l.isNumeric()) {
                EventActivityMap<T>::numeric_activity[l.variable()] += increment;
                need_rescaling = EventActivityMap<T>::numeric_activity[l.variable()] > maxActivity;
            } else {
                EventActivityMap<T>::boolean_activity[l.variable()] += increment;
                need_rescaling = EventActivityMap<T>::boolean_activity[l.variable()] > maxActivity;
                if(l.hasSemantic()) {
                    auto c{solver.boolean.getEdge(l)};
                    EventActivityMap<T>::numeric_activity[c.from] += increment;
                    EventActivityMap<T>::numeric_activity[c.to] += increment;
                    need_rescaling |= EventActivityMap<T>::numeric_activity[c.from] > maxActivity;
                    need_rescaling |= EventActivityMap<T>::numeric_activity[c.to] > maxActivity;
                }
            }
            return need_rescaling;
        }
            
        
        /**
         * Updates the activity map of the variables involved in the given clause, also applies the decay to each value
         * @tparam Clause
         * @param clause
         */
        template<class Iterable>
        void update(Iterable &clause, const Scheduler<T>& sched) noexcept {
          //            static_assert(traits::is_iterable_v<Clause>, "Expected
          //            iterable for Clause type");
          bool normalize = false;

          for (const auto gl : clause) {

#ifdef DEBUG_HEURISTICS
                std::cout << " increment activity because of " << sched.prettyLiteral(gl) << std::endl;
#endif
                
                if(LTYPE(gl) == BOUND_LIT) {
                    normalize |= incrementActivity(EVENT(sched.getBoundLiteral(FROM_GEN(gl))));
                } else {
                    DistanceConstraint<T> bc{sched.getEdge(FROM_GEN(gl))};
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

                auto [l, u] = std::ranges::minmax_element(this->numeric_activity);
                this->for_each([lb = *l, gap = *u - *l](auto &val) {
                  val = (val - lb) / gap * baseGap + baseIncrement;
                });
                increment = baseIncrement;
            }

            increment /= decay;
        }
        
        template<class Iterable>
        void update(Iterable &clause, const Solver<T>& solver) noexcept {

            bool normalize = false;

            
            for (const auto l : clause) {
                
//                if(solver.num_fails >= 20901)
//                    std::cout << " " << l ;
//                
#ifdef DEBUG_HEURISTICS
                std::cout << " increment activity because of " << l << std::endl;
#endif
                
                normalize |= incrementActivity(l, solver);
                
            }
            
//            if(solver.num_fails >= 20901)
//                std::cout << "\n--\n" ;
//            
#ifdef DEBUG_HEURISTICS
            for(event x{0}; x<static_cast<event>(EventActivityMap<T>::numeric_activity.size()); ++x) {
                std::cout << std::setw(6) << x;
            }
            std::cout << std::endl;
            for(event x{0}; x<static_cast<event>(EventActivityMap<T>::numeric_activity.size()); ++x) {
                std::cout << std::setw(6) << EventActivityMap<T>::numeric_activity[x];
            }
            std::cout << std::endl;
 #endif

            // protect against overflow
            if (normalize) {
                
//                if(solver.num_fails >= 20901)
//                std::cout << "\nnormalize (" << increment << ")\n" ;
                
#ifdef DEBUG_HEURISTICS
                std::cout << "\nnormalize (" << increment << ")\n" ;
#endif

                auto [nl, nu] = std::ranges::minmax_element(this->numeric_activity);
                auto [bl, bu] = std::ranges::minmax_element(this->boolean_activity);
                double l{std::numeric_limits<double>::max()};
                if(nl != this->numeric_activity.end())
                    l = *nl;
                if(bl != this->boolean_activity.end() and *bl < l)
                    l = *bl;
                
                double u{-std::numeric_limits<double>::max()};
                if(nu != this->numeric_activity.end())
                    u = *nu;
                if(bu != this->boolean_activity.end() and *bu > u)
                    u = *bu;
                
//                auto l{std::min(*nl,*bl)};
//                auto u{std::max(*nu,*bu)};
                this->for_each([lb = l, gap = u - l](auto &val) {
                  val = (val - lb) / gap * baseGap + baseIncrement;
                });
                
                increment = baseIncrement;
                
            }

            increment /= decay;
            
            
//            if(solver.num_fails >= 20901)
//                std::cout << "end update (" << solver.num_fails <<")\n";

            //            this->display(std::cout);
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

#endif // TEMPO_DECAYINGEVENTACTIVITYMAP_HPP

// x - y <= a
// x <= b
//
