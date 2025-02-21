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

#include "util/traits.hpp"

namespace detail {
    template<std::size_t Idx, typename Tuple>
    void printElem(std::ostream& os, const Tuple &tuple) {
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
auto printRange(R &&range, std::ostream &os) -> std::ostream& {
    os << "[";
    bool first = true;
    for (const auto &elem : std::forward<R>(range)) {
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



#endif //PRINTING_HPP
