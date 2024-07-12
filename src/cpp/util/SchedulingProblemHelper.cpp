/**
* @author Tim Luchterhand
* @date 01.07.24
* @brief
*/

#include "util/SchedulingProblemHelper.hpp"

namespace tempo{

    unsigned VarTaskMapping::operator()(var_t variable) const noexcept {
        return varToTask[variable - offset];
    }

    bool VarTaskMapping::contains(var_t variable) const noexcept {
        const var_t idx = variable - offset;
        return idx < varToTask.size() and varToTask[idx] != NoTask;
    }
}