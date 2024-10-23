/**
* @author Tim Luchterhand
* @date 23.10.24
* @brief
*/

#ifndef TEMPO_RELAXATION_POLICY_FACTORIES_HPP
#define TEMPO_RELAXATION_POLICY_FACTORIES_HPP
#include <ranges>

#include "util/factory_pattern.hpp"
#include "util/enum.hpp"
#include "util/traits.hpp"
#include "util/Options.hpp"
#include "Model.hpp"
#include "RelaxationPolicy.hpp"

namespace tempo::heuristics {

    PENUM(RelaxPolicy, RandomSubset, AllButOneResource, RandomResource)

    struct RelaxationPolicyParams {
        double relaxRatio;
        double ratioDecay;
    };

    template<concepts::scalar T, resource_expression R>
    using RelaxationPolicy = VariantPolicy<RandomSubset<T>,
            FixRandomDisjunctiveResource<R>, RelaxRandomDisjunctiveResource<R>>;

    MAKE_TEMPLATE_FACTORY(RandomSubset, ESCAPE(resource_range R),
                          ESCAPE(R &&resources, const RelaxationPolicyParams &params)) {
            auto variables = booleanVarsFromResources(std::forward<R>(resources));
            return RandomSubset(std::move(variables), params.relaxRatio, params.ratioDecay);
        }
    };

    MAKE_TEMPLATE_FACTORY(AllButOneResource, ESCAPE(resource_range R),
                          ESCAPE(R && resources, const RelaxationPolicyParams &)) {
            return FixRandomDisjunctiveResource(std::forward<R>(resources));
        }
    };

    MAKE_TEMPLATE_FACTORY(RandomResource, ESCAPE(resource_range R),
                          ESCAPE(R && resources, const RelaxationPolicyParams &)) {
            return RelaxRandomDisjunctiveResource(std::forward<R>(resources));
        }
    };

    template<concepts::scalar T, resource_expression R>
    MAKE_P_FACTORY_PATTERN(RelaxationPolciy, ESCAPE(RelaxationPolicy<T, R>),
                           RandomSubset, AllButOneResource, RandomResource)

    template<resource_range RR>
    auto make_relaxation_policy(RelaxPolicy type, RR &&resources, const RelaxationPolicyParams &params) {
        using R = std::ranges::range_value_t<RR>;
        using Factory = RelaxationPolciyFactory<tempo::detail::timing_type_from_resource_t<R>, R>;
        return Factory::getInstance().create(to_string(type), std::forward<RR>(resources), params);
    }
}

#endif //TEMPO_RELAXATION_POLICY_FACTORIES_HPP
