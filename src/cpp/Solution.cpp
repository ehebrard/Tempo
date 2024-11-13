/**
* @author Tim Luchterhand
* @date 12.11.24
* @brief
*/

#include "Solution.hpp"

namespace tempo {
    bool BooleanSolution::value(tempo::var_t x) const noexcept {
        return polarities[x];
    }

    BooleanSolution::BooleanSolution(std::vector<bool> polarities) noexcept: polarities(std::move(polarities)) {}

    std::size_t BooleanSolution::size() const noexcept {
        return polarities.size();
    }

    auto BooleanSolution::end() const noexcept -> std::vector<bool>::const_iterator {
        return polarities.cend();
    }

    auto BooleanSolution::begin() const noexcept -> std::vector<bool>::const_iterator {
        return polarities.cbegin();
    }
}