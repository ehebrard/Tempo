/**
* @author Tim Luchterhand
* @date 13.11.24
* @brief contains utilities for measuring distances between tasks (intervals)
*/

#ifndef TEMPO_TASK_DISTANCE_HPP
#define TEMPO_TASK_DISTANCE_HPP

#include <vector>
#include <concepts>

#include "util/traits.hpp"
#include "util/Matrix.hpp"
#include "edge_distance.hpp"
#include "DirectedGraph.hpp"
#include "Model.hpp"

namespace tempo {

    template<typename Arc>
    concept labeled_arc = concepts::same_template<Arc, LabeledEdge> or
                          concepts::same_template<Arc, StampedLabeledEdge>;

    namespace detail {
        template<typename T, typename Arc>
        void getDistances(Matrix<T> &matrix, const std::vector<std::vector<Arc>> &adjacency, bool incoming) noexcept {
            for (auto [src, neighbors] : iterators::enumerate(adjacency, 0)) {
                for (const auto &dest : neighbors) {
                    const auto s = incoming ? static_cast<int>(dest) : src;
                    const auto d = incoming ? src : static_cast<int>(dest);
                    matrix(s, d) = dest.label();
                }
            }
        }
    }

    /**
     * @brief Provides function like access to a bound provider for easy distance estimation between tasks
     * @tparam T timing type
     * @tparam BoundProvider boudn provider type
     */
    template<concepts::scalar T, bound_provider BoundProvider>
    class TaskDistanceView {
    protected:
        const std::vector<Interval<T>> &tasks;
        const BoundProvider &boundProvider;

        auto rawDistance(unsigned taskFrom, unsigned taskTo) const
                -> std::pair<decltype(std::declval<BoundProvider>().upper(0)), T> {
            if (taskFrom == taskTo) {
                return {0, 0};
            }

            const auto srcVar = tasks[taskFrom].start;
            const auto destVar = tasks[taskTo].end;
            return {boundProvider.upper(destVar.id()) - boundProvider.lower(srcVar.id()), destVar.offset()};
        }
    public:
        TaskDistanceView(const std::vector<Interval<T>> &tasks, const BoundProvider &boundProvider) noexcept:
            tasks(tasks), boundProvider(boundProvider) {}


        [[nodiscard]] auto operator()(unsigned taskFrom, unsigned taskTo) const {
            auto [dist, offset] = rawDistance(taskFrom, taskTo);
            return dist + offset;
        }
    };

    /**
     * @brief Provides a task distance view considering both information from a graph and from a bound provider.
     * @tparam T timing type
     * @tparam BoundProvider bound provider type
     */
    template<concepts::scalar T, bound_provider BoundProvider>
    class TaskDistanceFunction: private TaskDistanceView<T, BoundProvider> {
        Matrix<T> graphDistances;
    public:

        template<labeled_arc Arc>
        TaskDistanceFunction(const std::vector<Interval<T>> &tasks, const DirectedGraph<Arc> &varGraph,
                             const BoundProvider &boundProvider) :
             TaskDistanceView<T, BoundProvider>(tasks, boundProvider),
             graphDistances(varGraph.size(), varGraph.size(), InfiniteDistance<T>) {
            detail::getDistances(graphDistances, varGraph.forward(), false);
            detail::getDistances(graphDistances, varGraph.backward(), true);
        }

        [[nodiscard]] T operator()(unsigned taskFrom, unsigned taskTo) const {
            if (taskFrom == taskTo) {
                return 0;
            }

            const auto srcVar = this->tasks[taskFrom].start;
            const auto destVar = this->tasks[taskTo].end;
            const auto [boundDistance, _] = TaskDistanceView<T, BoundProvider>::rawDistance(taskFrom, taskTo);
            return std::min(graphDistances(srcVar.id(), destVar.id()), boundDistance) + destVar.offset();
        }
    };
}

#endif //TEMPO_TASK_DISTANCE_HPP
