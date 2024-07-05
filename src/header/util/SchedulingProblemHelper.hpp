/**
* @author Tim Luchterhand
* @date 01.07.24
* @brief utility for extracting timing information between tasks
*/

#ifndef TEMPO_SCHEDULINGPROBLEMHELPER_HPP
#define TEMPO_SCHEDULINGPROBLEMHELPER_HPP

#include <concepts>
#include <limits>
#include <Iterators.hpp>
#include <vector>
#include <ranges>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "DirectedGraph.hpp"
#include "Model.hpp"
#include "util/Matrix.hpp"
#include "util/traits.hpp"

namespace tempo {

    namespace detail {
        template<concepts::scalar T>
        struct NoValue {
            static constexpr auto value() {
                if constexpr (std::floating_point<T>) {
                    return std::numeric_limits<T>::infinity();
                } else {
                    return std::numeric_limits<T>::max();
                }
            }
        };

        template<concepts::scalar T>
        inline constexpr auto NoDistance = detail::NoValue<T>::value();

        template<typename T, typename Arc>
        void getDistances(Matrix<T> &matrix, const std::vector<std::vector<Arc>> &adjacency) noexcept {
            for (auto [src, neighbors] : iterators::enumerate(adjacency)) {
                for (const auto &dest : neighbors) {
                    matrix(src, dest) = dest.label();
                }
            }
        }

        template<typename Arc>
        concept labeled_arc = concepts::same_template<Arc, LabeledEdge> or
                              concepts::same_template<Arc, StampedLabeledEdge>;

        template<typename P, typename T>
        concept bound_provider = requires(const P instance, var_t var) {
            { instance.upper(var) } -> std::same_as<T>;
            { instance.lower(var) } -> std::same_as<T>;
        };

        template<concepts::scalar T, bound_provider<T> BoundProvider>
        class TaskDistanceFunction {
            const std::vector<Interval<T>> &tasks;
            const BoundProvider &boundProvider;
            Matrix<T> graphDistances;
        public:

            template<labeled_arc Arc>
            TaskDistanceFunction(const std::vector<Interval<T>> &tasks, const DirectedGraph<Arc> &varGraph,
                                 const BoundProvider &boundProvider) :
                    tasks(tasks), boundProvider(boundProvider),
                    graphDistances(varGraph.vertexCount(), varGraph.vertexCount(), NoDistance<T>) {
                getDistances(graphDistances, varGraph.forward());
                getDistances(graphDistances, varGraph.backward());
            }

            [[nodiscard]] T operator()(unsigned taskFrom, unsigned taskTo) const {
                const auto srcVar = tasks[taskFrom].start;
                const auto destVar = tasks[taskTo].end;
                const auto boundDistance = boundProvider.upper(srcVar) - boundProvider.lower(destVar);
                return std::min(graphDistances(srcVar, destVar), boundDistance) + destVar.offset();
            }
        };
    }

    /**
     * Mapping from variable id to task id
     */
    class VarTaskMapping {
        std::vector<unsigned> varToTask;
        var_t offset;
        static constexpr unsigned NoTask = std::numeric_limits<unsigned>::max();
    public:
        /**
         * Ctor
         * @tparam Tasks range containing tempo::Interval
         * @param tasks all tasks of the problem
         */
        template<concepts::ttyped_range<tempo::Interval> Tasks>
        requires(std::ranges::sized_range<Tasks>)
        constexpr explicit VarTaskMapping(const Tasks &tasks) : varToTask(2 * std::ranges::size(tasks), NoTask),
                                                                offset(std::numeric_limits<var_t>::max()) {

            var_t max = 0;
            for (const auto &task : tasks) {
                max = std::max({max, task.start.id(), task.end.id()});
                offset = std::min({offset, task.start.id(), task.end.id()});
            }

            if (not varToTask.empty() and max - offset >= 2 * std::ranges::size(tasks)) {
                throw std::runtime_error("expected continuous ranges of variable ids");
            }

            for (auto [idx, task]: iterators::enumerate(tasks)) {
                varToTask[task.start.id() - offset] = idx;
                varToTask[task.end.id() - offset] = idx;
            }
        }

        /**
         * Mapping from variable id to task
         * @param variable id of the variable
         * @return corresponding task id
         */
        [[nodiscard]] unsigned operator()(var_t variable) const noexcept;

        /**
         * Whether the mapping contains a variable
         * @param variable id of the variable
         * @return true if variable is in mapping, false otherwise
         */
        [[nodiscard]] bool contains(var_t variable) const noexcept;

        /**
         * @return number of tasks in the mapping
         */
        [[nodiscard]] std::size_t size() const noexcept;

    };

    template<typename T>
    class Solver;

    /**
     * Class that represents a scheduling problem and provides helper functions
     */
    template<concepts::scalar T, SchedulingResource R>
    class SchedulingProblemHelper {
        using TaskVec = std::vector<Interval<T>>;
        using RVec = std::vector<R>;
        using PrecVec = std::vector<DistanceConstraint<T>>;
        TaskVec t;
        RVec res;
        PrecVec precs;
        VarTaskMapping v2t;

    public:
        /**
         * CTor
         * @param tasks tasks contained in the problem
         * @param resources resources in the problem
         * @param precedences precedences between variables in the problem
         */
        SchedulingProblemHelper(TaskVec tasks, RVec resources, PrecVec precedences) :
                t(std::move(tasks)), res(std::move(resources)), precs(std::move(precedences)), v2t(t) {}

        /**
         * read only access to the tasks
         * @return vector with the tasks
         */
        [[nodiscard]] auto tasks() const noexcept -> const TaskVec & {
            return t;
        }

        /**
         * read only access to the resources
         * @return vector with the resources
         */
        [[nodiscard]] auto resources() const noexcept -> const RVec & {
            return res;
        }

        /**
         * read only access to the precedences
         * @return vector with the precedences
         */
        [[nodiscard]] auto precedences() const noexcept -> const PrecVec & {
            return precs;
        }

        /**
         * @return task associated with the variable
         * @param variable id of a variable
         * @note UB if variable not associated with a task
         */
        [[nodiscard]] auto getTask(var_t variable) const noexcept -> TaskVec::const_reference {
            return t[v2t(variable)];
        }

        /**
         * Whether a given task id belongs to the problem
         * @param taskId id of a task
         * @return true if task id exists in problem, false otherwise
         */
        [[nodiscard]] bool hasTask(unsigned taskId) const noexcept {
            return taskId < t.size();
        }

        /**
         * Whether a given variable id is associated with a task
         * @param variable id of a variable
         * @return whether the variable is contained in the variable-to-task mapping
         */
        [[nodiscard]] bool hasVariable(var_t variable) const noexcept {
            return v2t.contains(variable);
        }

        /**
         * Gets a functor object that represents temporal distances between tasks. The distance between two tasks t1 and
         * t2 is equal to dist(start_t1, end_t2)
         * @tparam BoundProvider Class that provides upper and lower bounds for numeric variables
         * @tparam Arc Arc type that contains a distance label
         * @param varGraph directed graph representing distances between variables
         * @param boundProvider object that provides upper and lower bounds
         * @return distance function object that provides distances between tasks
         */
        template<detail::bound_provider<T> BoundProvider, detail::labeled_arc Arc>
        [[nodiscard]] auto getTaskDistances(const DirectedGraph<Arc> &varGraph,
                                            const BoundProvider &boundProvider) const {
            return detail::TaskDistanceFunction<T, BoundProvider>(t, varGraph, boundProvider);
        }

        /**
         * Function overload of getTaskDistances
         * @param solver solver instance used to construct the task distance function
         * @return distance function object that provides distances between tasks
         */
        [[nodiscard]] auto getTaskDistances(const Solver<T> &solver) const {
            using BP = std::remove_cvref_t<decltype(solver.numeric)>;
            return detail::TaskDistanceFunction<T, BP>(t, solver.core, solver.numeric);
        }
    };
}

#endif //TEMPO_SCHEDULINGPROBLEMHELPER_HPP
