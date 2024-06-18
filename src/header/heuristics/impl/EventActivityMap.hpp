//
// Created by tluchterha on 22/11/22.
//

#ifndef SCHEDCL_EVENTACTIVITYMAP_HPP
#define SCHEDCL_EVENTACTIVITYMAP_HPP

#include <concepts>

//#include "Scheduler.hpp"

namespace tempo {
 
template<typename T>
class Solver;
}

namespace tempo::heuristics::impl {
    /**
     * @brief Class that can be used to record the activity on distance constraints
     */
template<typename T>
    class EventActivityMap {
    public:
        /**
         * CTor. Initializes activity of all literals with 1.
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to construct the ActivityMap
         */
        
        explicit EventActivityMap(Solver<T> &solver)
//        : sched(scheduler)
        {
//            numNodes = scheduler.numEvent();
            numeric_activity.resize(solver.numeric.size(), 1);
            boolean_activity.resize(solver.boolean.size(), 1);
            solver.setActivityMap(this);
        }

        
        
        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param edge constraint
         * @return
         */
        constexpr double get(const DistanceConstraint<T>& c) const noexcept {
            return numeric_activity[c.from] + numeric_activity[c.to];
        }
        
        constexpr double get(const Literal<T> l,
                             const Solver<T> &solver) const noexcept {
          double a{0};

          if (l.isNumeric()) {
            a = numeric_activity[l.variable()];
          } else {
            a = boolean_activity[l.variable()];
            if (l.hasSemantic()) {
              a += get(solver.boolean.getEdge(l));
              a += get(solver.boolean.getEdge(~l));
            }
          }
          return a;
        }

        constexpr double get(const var_t x, const Solver<T>& solver) const noexcept {
          double a{boolean_activity[x]};
          //            double a{0};
          if (solver.boolean.hasSemantic(x)) {
              a += get(solver.boolean.getEdge(true, x));// / 1000;
              a += get(solver.boolean.getEdge(false, x));// / 1000;
          }

            return a;
        }

//        /**
//         * Gets the activity for a given choice point
//         * @tparam T
//         * @param choicePoint
//         * @return
//         */
//        constexpr double &get(const var x) noexcept {
//            DistanceConstraint<T> left{sched.getEdge(POS(x))};
//            DistanceConstraint<T> right{sched.getEdge(NEG(x))};
//            return activity[left.from] + activity[left.to] + activity[right.from] + activity[right.to];
//        }

        /**
         * Applies the given functor to all entries in the activity map
         * @param functor
         */
        void for_each(const std::invocable<double &> auto &functor) {
            for (auto &val : numeric_activity) {
                functor(val);
            }
            for (auto &val : boolean_activity) {
                functor(val);
            }
        }
        
        /**
         * Displays weights
         */
//        template<typename T>
        std::ostream& display(std::ostream& os) {
          //            for(auto i{0}; i<numeric_activity.size(); ++i) {
          //                if(numeric_activity[i] > 1)
          //                    os << " x" << i << ":" << numeric_activity[i];
          //            }
          //            for(auto i{0}; i<boolean_activity.size(); ++i) {
          //                if(boolean_activity[i] > 1)
          //                    os << " b" << i << ":" << boolean_activity[i];
          //            }
          for (size_t i{0}; i < numeric_activity.size(); ++i) {
            os << std::setprecision(7) << std::setw(10) << numeric_activity[i];
          }
          os << " |";
          for (size_t i{0}; i < boolean_activity.size(); ++i) {
            os << std::setprecision(7) << std::setw(10) << boolean_activity[i];
          }
          os << "\n";
          return os;
        }

    protected:

//        Scheduler<T>& sched;
        std::vector<double> numeric_activity{};
        std::vector<double> boolean_activity{};
    };
}

#endif //SCHEDCL_EVENTACTIVITYMAP_HPP
