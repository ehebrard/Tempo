/**
* @author Tim Luchterhand
* @date 23.10.24
* @brief Contains factory methods for different standard relaxation policies defined in heuristics/RelaxationPolicy.hpp
*/

#ifndef TEMPO_RELAXATION_POLICY_FACTORIES_HPP
#define TEMPO_RELAXATION_POLICY_FACTORIES_HPP
#include <ranges>

#include "util/factory_pattern.hpp"
#include "util/enum.hpp"
#include "util/traits.hpp"
#include "util/Options.hpp"
#include "Model.hpp"
#include "relaxation_policies.hpp"

namespace tempo::lns {

    /**
     * @brief Relaxation policy type
     */
    PENUM(RelaxPolicy, RandomSubset, AllButOneResource, RandomResource, RandomTasks, Chronologically)

    /**
     * @brief Configuration parameters for relaxation policies
     */
    struct RelaxationPolicyParams {
        double relaxRatio; ///< percentage of variables to relax (RandomSubset policy)
        double ratioDecay; ///< relaxation ratio decay
        unsigned numScheduleSlices; ///< number of schedule slices for chronological task relaxation
    };

    // --- add further relaxation policies to this type ---
    template<concepts::scalar T, resource_expression R>
    using RelaxationPolicy = VariantPolicy<RandomSubset<T>, FixRandomDisjunctiveResource<R>,
            RelaxRandomDisjunctiveResource<R>, RelaxTasks<T>, RelaxChronologically<T>>;

    // --- Define policy factories here

    MAKE_FACTORY(RandomSubset, ESCAPE(auto &&, resource_range auto &&resources, const RelaxationPolicyParams &params)) {
            auto variables = booleanVarsFromResources(std::forward<decltype(resources)>(resources));
            return RandomSubset(std::move(variables), params.relaxRatio, params.ratioDecay);
        }
    };

    MAKE_FACTORY(AllButOneResource, ESCAPE(auto &&, resource_range auto &&resources, const RelaxationPolicyParams &)) {
            return FixRandomDisjunctiveResource(std::forward<decltype(resources)>(resources));
        }
    };

    MAKE_FACTORY(RandomResource, ESCAPE(auto &&, resource_range auto && resources, const RelaxationPolicyParams &)) {
            return RelaxRandomDisjunctiveResource(std::forward<decltype(resources)>(resources));
        }
    };

    MAKE_TEMPLATE_FACTORY(RandomTasks, ESCAPE(concepts::scalar T, resource_range R),
                          ESCAPE(std::vector<Interval<T>> tasks, R &&resources, const RelaxationPolicyParams &params)) {
            return RelaxTasks(std::move(tasks), std::forward<R>(resources), params.relaxRatio, params.ratioDecay);
        }
    };

    MAKE_TEMPLATE_FACTORY(Chronologically, ESCAPE(concepts::scalar T, resource_range R),
                          ESCAPE(std::vector<Interval<T>> tasks, R &&resources, const RelaxationPolicyParams &params)) {
            return RelaxChronologically(std::move(tasks), std::forward<R>(resources), params.numScheduleSlices);
        }
    };

    // --- Don't forget to add the new type here too (at the end of the variadic argument list)

    template<concepts::scalar T, resource_expression R>
    MAKE_P_FACTORY_PATTERN(RelaxationPolciy, ESCAPE(RelaxationPolicy<T, R>),
                           RandomSubset, AllButOneResource, RandomResource, RandomTasks, Chronologically)

    /**
     * Factory method for relaxation policies
     * @tparam RR type of resource range
     * @param type relaxation policy type
     * @param resources resource expressions of the problem
     * @param params policy parameters
     * @return relaxation policy for large neighborhood search
     */
    template<concepts::scalar T, resource_range RR>
    auto make_relaxation_policy(RelaxPolicy type, std::vector<Interval<T>> tasks, RR &&resources,
                                const RelaxationPolicyParams &params) {
        using R = std::ranges::range_value_t<RR>;
        using Factory = RelaxationPolciyFactory<tempo::detail::timing_type_from_resource_t<R>, R>;
        return Factory::getInstance().create(to_string(type), std::move(tasks), std::forward<RR>(resources), params);
    }
}

#endif //TEMPO_RELAXATION_POLICY_FACTORIES_HPP
