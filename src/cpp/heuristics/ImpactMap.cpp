/**
* @author Tim Luchterhand
* @date 11.03.25
* @file ImpactMap.cpp
* @brief
*/

#include "heuristics/impl/ImpactMap.hpp"

namespace tempo::heuristics::impl {
    void ImpactMap::updateImpact(unsigned numFixed) {
        auto &impactVal = impact.at(litId);
        impactVal = (impactVal + static_cast<ImpactMeasure>(numFixed)) / 2;
    }

    auto ImpactMap::get(info_t lit) const -> ImpactMeasure {
        return impact.at(lit);
    }

    auto ImpactMap::getMap() const noexcept -> const std::vector<ImpactMeasure> & {
        return impact;
    }
}
