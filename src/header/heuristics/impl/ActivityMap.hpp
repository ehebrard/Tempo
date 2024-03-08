//
// Created by tluchterha on 22/11/22.
//

#ifndef SCHEDCL_ACTIVITYMAP_HPP
#define SCHEDCL_ACTIVITYMAP_HPP

#include <concepts>

#include "DistanceMatrix.hpp"

namespace schedcl {
    template<typename T>
    class Scheduler;
}

namespace schedcl::heuristics::impl {
    /**
     * @brief Class that can be used to record the activity on distance constraints
     */
    class ActivityMap {
    public:
        /**
         * CTor. Initializes activity of all literals with 1.
         * @tparam T type of scheduler
         * @param scheduler scheduler for which to construct the ActivityMap
         */
        template<typename T>
        explicit ActivityMap(const Scheduler<T> &scheduler) {
            using CPSequence = std::remove_cvref_t<decltype(scheduler.getChoicePoints())>;
            static_assert(traits::is_iterable_v<CPSequence>);
            static_assert(traits::is_same_template_v<DistanceConstraint, traits::value_type_t<CPSequence>>);

            event minEvent = std::numeric_limits<event>::max();
            event maxEvent = std::numeric_limits<event>::lowest();
            for (const auto &cp : scheduler.getChoicePoints()) {
                minEvent = std::min({minEvent, cp.from, cp.to});
                maxEvent = std::max({minEvent, cp.from, cp.to});
            }

            assert(minEvent >= 0);
            indexOffset = minEvent;
            numNodes = maxEvent - minEvent + 1;
            activity.resize(numNodes * (numNodes - 1) / 2, 1);
        }

        /**
         * Checks whether a given constraint is registered in the activity map
         * @tparam T type of constraint
         * @param constraint constraint to check
         * @return true if activity map contains a value for the given constraint, false otherwise
         */
        template<typename T>
        constexpr bool contains(const DistanceConstraint<T> &constraint) const noexcept {
            const auto lower = std::min(constraint.from, constraint.to);
            const auto upper = std::max(constraint.from, constraint.to);
            return lower >= 0 && upper >= 0 && static_cast<std::size_t>(lower) >= indexOffset &&
                   static_cast<std::size_t>(upper) < numNodes + indexOffset;
        }

        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param choicePoint
         * @return
         */
        template<typename T>
        constexpr double get(const DistanceConstraint<T> & choicePoint) const noexcept {
            assert(contains(choicePoint));
            return activity[getIndex(choicePoint)];
        }

        /**
         * Gets the activity for a given choice point
         * @tparam T
         * @param choicePoint
         * @return
         */
        template<typename T>
        constexpr double &get(const DistanceConstraint<T> & choicePoint) noexcept {
            assert(contains(choicePoint));
            return activity[getIndex(choicePoint)];
        }

        /**
         * Applies the given functor to all entries in the activity map
         * @param functor
         */
        void for_each(const std::invocable<double &> auto &functor) {
            for (auto &val : activity) {
                functor(val);
            }
        }

        /**
         * Displays weights
         */
        template <typename T> std::ostream &display(std::ostream &os) {
          os << "     ";
          for (auto i{0}; i < numNodes; ++i) {
            os << std::setw(5) << i;
          }
          for (auto i{0}; i < numNodes; ++i) {
            os << std::setw(3) << i << ": ";
            for (auto j{0}; j < numNodes; ++j)
              if (i != j) {
                DistanceConstraint<T> x{i, j, 0};
                os << std::setw(5) << getIndex(x);
                //                    activity[getIndex(x)];
              } else {
                os << std::setw(5) << " ";
              }
            os << std::endl;
          }
          return os;
        }

    protected:
        /**
         * Returns the index to the activity map for a given choice point
         * @tparam T
         * @param choicePoint
         * @return
         * @note Elements are stored in an upper triangular matrix with row major organization. Also, direction of
         * distance constraints is ignored. This means, that for example [4 -> 6] has the same index as [6 -> 4]
         */
        template<typename T>
        constexpr std::size_t getIndex(const DistanceConstraint<T> &choicePoint) const noexcept {

            const auto row = std::min(choicePoint.from, choicePoint.to) - indexOffset;
            const auto col = std::max(choicePoint.from, choicePoint.to) - indexOffset;
            const auto numElementsInRows = row * numNodes - row * (row + 1) / 2;

            //            std::cout << std::endl << "get index " <<
            //            choicePoint.from << " / "
            //            << choicePoint.to << " row=" << row << " col=" << col
            //            << " numElementsInRows=" << numElementsInRows
            //            << " => " << (numElementsInRows + col - (row + 1)) <<
            //            std::endl;

            return numElementsInRows + col - (row + 1);
        }

        std::vector<double> activity{};
        std::size_t numNodes{};
        std::size_t indexOffset{};
    };
}

#endif //SCHEDCL_ACTIVITYMAP_HPP
