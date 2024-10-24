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
#include "util/distance.hpp"
#include "Literal.hpp"

namespace tempo {

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
                    graphDistances(varGraph.size(), varGraph.size(), InfiniteDistance<T>) {
                getDistances(graphDistances, varGraph.forward(), false);
                getDistances(graphDistances, varGraph.backward(), true);
            }

            [[nodiscard]] T operator()(unsigned taskFrom, unsigned taskTo) const {
                if (taskFrom == taskTo) {
                    return 0;
                }

                const auto srcVar = tasks[taskFrom].start;
                const auto destVar = tasks[taskTo].end;
                const auto boundDistance = boundProvider.upper(destVar.id()) - boundProvider.lower(srcVar.id());
                return std::min(graphDistances(srcVar.id(), destVar.id()), boundDistance) + destVar.offset();
            }
        };
    }

    template<typename T>
    concept SchedulingResource = requires(const T instance, unsigned taskId) {
        { instance.tasks() } -> concepts::typed_range<unsigned>;
        { instance.tasks() } -> std::ranges::sized_range;
        { instance.resourceCapacity() } -> concepts::scalar;
        { instance.getDemand(taskId) } -> concepts::scalar;
    };

    /**
     * Mapping from variable id to task id
     */
    class VarTaskMapping {
        std::vector<unsigned> varToTask;
        var_t offset{};
        static constexpr unsigned NoTask = std::numeric_limits<unsigned>::max();
    public:
        static constexpr unsigned NumTemporalVarPerTask = 3;
        VarTaskMapping() = default;

        /**
         * Ctor
         * @tparam Tasks range containing tempo::Interval
         * @param tasks all tasks of the problem
         */
        template<concepts::ttyped_range<tempo::Interval> Tasks>
        requires(std::ranges::sized_range<Tasks>)
        constexpr explicit VarTaskMapping(const Tasks &tasks) :
                varToTask(NumTemporalVarPerTask * std::ranges::size(tasks), NoTask),
                offset(std::numeric_limits<var_t>::max()) {

            var_t max = 0;
            for (const auto &task : tasks) {
                max = std::max({max, task.start.id(), task.end.id()});
                offset = std::min({offset, task.start.id(), task.end.id()});
            }

            if (not varToTask.empty() and max - offset >= NumTemporalVarPerTask * std::ranges::size(tasks)) {
                throw std::runtime_error("expected continuous ranges of variable ids");
            }

            for (auto [idx, task]: iterators::enumerate(tasks)) {
                auto &start = varToTask[task.start.id() - offset];
                auto &end = varToTask[task.end.id() - offset];
                if (start != NoTask or end != NoTask) {
                    throw std::runtime_error("multiple tasks with same variables");
                }

                start = idx;
                end = idx;
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
        Interval<T> sched;

    public:
        using Time = T;
        using Resource = R;
        SchedulingProblemHelper() = default;

        /**
         * CTor
         * @param tasks tasks contained in the problem
         * @param resources resources in the problem
         * @param precedences precedences between variables in the problem
         * @param schedule Interval that represents the entire schedule
         */
        SchedulingProblemHelper(TaskVec tasks, RVec resources, PrecVec precedences, Interval<T> schedule) :
                t(std::move(tasks)), res(std::move(resources)), precs(std::move(precedences)), v2t(t), sched(schedule) {}

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
         * gets the mapping function from tasks to variables
         * @return mapping from variables to tasks
         */
        [[nodiscard]] auto getMapping() const noexcept -> VarTaskMapping {
            return v2t;
        }

        /**
         * Access to the interval representing the entire schedule
         */
        auto schedule() const noexcept -> const Interval<T> & {
            return sched;
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

        /**
         * Gets all non-fixed literals that have a semantic and correspond to task edges
         * @param solver solver from which to extract the literals
         * @return list of free literals corresponding to edges between tasks
         */
        auto getSearchLiterals(const Solver<T> &solver) const -> std::vector<Literal<T>> {
            std::vector<Literal<T>> literals;
            for (auto var : solver.getBranch()) {
                if (solver.boolean.hasSemantic(var)) {
                    auto edge = solver.boolean.getEdge(true, var);
                    if (hasVariable(edge.from) and hasVariable(edge.to)) {
                        literals.emplace_back(solver.boolean.getLiteral(true, var));
                    }
                }
            }

            return literals;
        }
    };
}

#endif //TEMPO_SCHEDULINGPROBLEMHELPER_HPP
