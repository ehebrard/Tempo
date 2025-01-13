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

template<std::ranges::range R>
auto printRange(const R &range, std::ostream &os) -> std::ostream& {
    os << "[";
    bool first = true;
    for (const auto &elem : range) {
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
