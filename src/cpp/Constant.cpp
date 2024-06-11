
#include "Constant.hpp"

using namespace tempo;

Explanation Constant::NoReason = Explanation(new Explainer(), NoHint);

const index_t Constant::NoIndex = static_cast<index_t>(-1);

const hint Constant::FactHint = static_cast<index_t>(-1);

const hint Constant::DecisionHint = static_cast<index_t>(0);

const info_t Constant::NoSemantic = static_cast<info_t>(0);

const var_t Constant::NoVarx = static_cast<var_t>(-1);

const index_t Constant::InfIndex = 0;

//const index_t Constant::IndexOfMax = 1;
