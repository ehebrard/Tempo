/**
* @author Tim Luchterhand
* @date 09.01.25
* @file printing.hpp
* @brief
*/

#ifndef PRINTING_HPP
#define PRINTING_HPP

#include <ranges>
#include <ostream>
#include <tuple>
#include <utility>
#include <sstream>
#include <vector>
#include <Iterators.hpp>

#include "util/traits.hpp"
#include "Solver.hpp"
#include "Model.hpp"

namespace detail {
    template<std::size_t Idx, typename Tuple>
    void printElem(std::ostream &os, const Tuple &tuple) {
        if constexpr (Idx != 0) {
            os << ", ";
        }

        os << std::get<Idx>(tuple);
    }
}

namespace std {
    template<tempo::concepts::printable ...Ts>
    std::ostream &operator<<(std::ostream &os, const std::tuple<Ts...> &tuple) {
        os << "(";
        [&]<std::size_t ... Idx>(std::index_sequence<Idx...>) {
            (detail::printElem<Idx>(os, tuple), ...);
        }(std::make_index_sequence<sizeof...(Ts)>());
        os << ")";
        return os;
    }

    template<tempo::concepts::printable T, tempo::concepts::printable U>
    std::ostream &operator<<(std::ostream &os, const std::pair<T, U> &pair) {
        os << "(" << pair.first << ", " << pair.second << ")";
        return os;
    }
}

template<std::ranges::range R> requires(tempo::concepts::printable<std::ranges::range_value_t<R>>)
auto printRange(R &&range, std::ostream &os) -> std::ostream & {
    os << "[";
    bool first = true;
    for (const auto &elem: std::forward<R>(range)) {
        if (first) {
            first = false;
        } else {
            os << ", ";
        }
        os << elem;
    }

    os << "]";
    return os;
}

template<typename T>
std::string prettyJob(const tempo::Interval<T> &task, const tempo::Solver<T> &S, const bool dur_flag) {
    std::stringstream ss;
    if (not S.numeric.hasSolution()) {
        throw std::runtime_error("Cannot print jobs, no solution found");
    }

    auto est{S.numeric.solutionLower(task.start)};
    auto lst{S.numeric.solutionUpper(task.start)};
    auto ect{S.numeric.solutionLower(task.end)};
    auto lct{S.numeric.solutionUpper(task.end)};
    if (S.boolean.value(task.exist)) {
        ss << "[";
        if (est == lst) {
            ss << est;
        } else {
            ss << est << "-" << lst;
        }
        ss << "..";
        if (ect == lct) {
            ss << ect;
        } else {
            ss << ect << "-" << lct;
        }
        ss << "]";

        if (dur_flag) {
            auto pmin{S.numeric.solutionLower(task.duration)};
            auto pmax{S.numeric.solutionUpper(task.duration)};
            ss << " (" << pmin << "-" << pmax << ")";
        }
    } else {
        ss << "removed";
    }

    return ss.str();
}


// TODO use a better resource interface (e.g. the one in parsing/scheduling_collection.hpp)
template<typename T>
void printResources(const tempo::Solver<T> &S, const std::vector<tempo::Interval<T>> &intervals,
                    const std::vector<std::vector<size_t>> &resource_tasks,
                    const std::vector<std::vector<std::vector<T>>> &resource_transitions = {}) {
    using namespace tempo;
    if (not S.numeric.hasSolution()) {
        throw std::runtime_error("Cannot print jobs, no solution found");
    }
    int i{0};
    for (auto &tasks: resource_tasks) {
        if (tasks.size() > 1) {
            std::vector<index_t> order;
            for (index_t j{0}; j < static_cast<index_t>(tasks.size()); ++j) {
                if (S.boolean.value(intervals[tasks[j]].exist)) {
                    order.push_back(j);
                }
            }

            std::sort(order.begin(), order.end(),
                      [&](const index_t a, const index_t b) {
                          return S.numeric.solutionLower(intervals[tasks[a]].start) < S.numeric.solutionLower(
                                     intervals[tasks[b]].start);
                      });

            std::cout << "resource " << i << ":";
            index_t n{Constant::NoIndex};
            for (auto o: order) {
                if (not resource_transitions.empty() and n != Constant::NoIndex) {
                    std::cout << " -> " << resource_transitions[i][order[n]][o] << " -> ";
                    ++n;
                } else {
                    n = 0;
                }

                std::cout << " job" << o << ":" << prettyJob(intervals[tasks[o]], S, false);
            }
            std::cout << std::endl;
        }

        ++i;
    }
}

template<tempo::concepts::scalar T, tempo::concepts::typed_range<tempo::Interval<T>> Tasks>
void printJobs(const tempo::Solver<T> &S, Tasks &&intervals) {
    for (auto [i, task]: iterators::enumerate(std::forward<Tasks>(intervals))) {
        std::cout << "job" << ++i << " (" << task.id() << "): " << prettyJob(task, S, true) << std::endl;
    }
}


#endif //PRINTING_HPP
