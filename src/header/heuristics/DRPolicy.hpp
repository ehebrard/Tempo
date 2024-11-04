/**
* @author Tim Luchterhand
* @date 10.10.24
* @brief Destroy-Release Policy interface and implementation
*/

#ifndef TEMPO_DRPOLICY_HPP
#define TEMPO_DRPOLICY_HPP

#include <concepts>
#include <vector>

#include "Literal.hpp"
#include "util/traits.hpp"
#include "RelaxationInterface.hpp"
#include "RelaxationPolicy.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {
    /**
     * @brief Destroy policy interface
     */
    template<typename T>
    concept destroy_policy = requires(T instance, unsigned numFail, AssumptionProxy<int> &proxy) {
        instance.notifySuccess(numFail);
        instance.notifyFailure(numFail);
        instance.requestNewRegion();
        { instance.relax(proxy) } -> concepts::ttyped_range<Literal>;
        { instance.newRegion() } -> std::convertible_to<bool>;
    };

    /**
     * @brief Repair policy interface
     */
    template<typename T>
    concept repair_policy = requires(T instance, unsigned numFail, AssumptionProxy<int> &proxy,
                                     std::vector<Literal<int>> region) {
        instance.notifySuccess(numFail);
        instance.notifyFailure(numFail);
        instance.notifyNewRegion(region);
        instance.fix(proxy);
        { instance.exhausted() } -> std::convertible_to<bool>;
    };

    /**
     * @brief Generic Repair-Distroy policy implementation
     * @tparam D destroy policy
     * @tparam R repair policy
     */
    template<destroy_policy D, repair_policy R>
    class RDPolicy {
        static constexpr auto Resolution = 100000;
        D destroy;
        R repair;
        double exhaustionProbability = 0;
    public:

        /**
         * Ctor
         * @tparam ArgD destroy policy type
         * @tparam ArgR repair policy type
         * @param destroyer destroy policy
         * @param repairer repair policy
         * @param exhaustionProbability probability of choosing a new region even if repair policy is not exhausted
         */
        template<destroy_policy ArgD, repair_policy ArgR>
        RDPolicy(ArgD &&destroyer, ArgR &&repairer, double exhaustionProbability = 0):destroy(
                std::forward<ArgD>(destroyer)), repair(std::forward<ArgR>(repairer)), exhaustionProbability(
                exhaustionProbability) {
            if (exhaustionProbability < 0 or exhaustionProbability > 1) {
                throw std::runtime_error("invalid exhaustion probability");
            }
        }

        /**
         * notifies both policy of a successful search
         * @param numFails number of fails of the solver after the search
         */
        void notifySuccess(unsigned numFails) {
            destroy.notifySuccess(numFails);
            repair.notifySuccess(numFails);
        }

        /**
         * notifies both policy of an  unsuccessful search
         * @param numFails number of fails of the solver after the search
         */
        void notifyFailure(unsigned numFails) {
            destroy.notifyFailure(numFails);
            repair.notifyFailure(numFails);
        }

        /**
         * relaxation policy interface
         * @tparam AI assumption interface
         * @param s solver proxy
         */
        template<assumption_interface AI>
        void relax(AI &s) {
            bool randomlyExhausted = randomEventOccurred<Resolution>(exhaustionProbability);
            if (repair.exhausted() or randomlyExhausted) {
                if (s.getSolver().getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- " << (randomlyExhausted ? "randomly " : "");
                    std::cout << "choosing new region" << std::endl;
                }

                destroy.requestNewRegion();
            }

            decltype(auto) fixedRegion = destroy.relax(s);
            if (s.getState() == AssumptionState::Fail) {
                return;
            }

            if (destroy.newRegion()) {
                repair.notifyNewRegion(std::forward<decltype(fixedRegion)>(fixedRegion));
            }

            repair.fix(s);
        }
    };

    /**
     * Helper method for constructing DR-policies facilitating template argument deduction
     * @tparam D destroy policy type
     * @tparam R repair policy type
     * @param destroy destroy policy
     * @param repair repair policy
     * @param exhaustionProbability probability of choosing a new region even if repair policy is not exhausted
     * @return constructed DR-policy
     */
    template<destroy_policy D, repair_policy R>
    auto make_RD_policy(D &&destroy, R &&repair, double exhaustionProbability = 0) {
        return RDPolicy<D, R>(std::forward<D>(destroy), std::forward<R>(repair), exhaustionProbability);
    }

    /**
     * @brief Generic destroy policy wrapper that transforms any relaxation policy into a destroy policy.
     * @details @copybrief
     * Caches and repeats the assumptions of the base relaxation policy and until a new region is requested
     * @tparam T timing type
     * @tparam R base relaxation policy used to make assumptions
     */
    template<concepts::scalar T, relaxation_policy R>
    class GenericDestroyPolicy : protected R {
        std::vector<Literal<T>> assumptionCache;
        bool newRegionRequested = true;
        bool isNewRegion = true;
    public:
        using R::notifyFailure;
        using R::notifySuccess;

        /**
         * Ctor
         * @tparam Args argument types
         * @param args arguments to base policy
         */
        template<typename ...Args>
        GenericDestroyPolicy(Args &&...args): R(std::forward<Args>(args)...) {}

        /**
         * Destroy policy interface. Calls the base policy to select a new region if a new region was requested.
         * Otherwise uses the cached assumptions to repeatedly explore the same region
         * @tparam AI assumption interface
         * @param proxy solver proxy
         */
        template<assumption_interface AI>
        auto relax(AI &proxy) -> const std::vector<Literal<T>> & {
            if (newRegionRequested) {
                newRegionRequested = false;
                AssumptionCollector<T, AI> collector(proxy);
                R::relax(collector);
                assumptionCache = std::move(collector.getAssumptions());
                isNewRegion = true;
            } else {
                proxy.makeAssumptions(assumptionCache);
                isNewRegion = false;
            }

            return assumptionCache;
        }

        /**
         * requests the policy to choose a new region
         */
        void requestNewRegion() {
            newRegionRequested = true;
        }

        /**
         * Indicates whether a new region has been chosen
         * @return true if a new region was chosen in the last call to relax, false otherwise
         */
        [[nodiscard]] bool newRegion() const noexcept {
            return isNewRegion;
        }
    };
}

#endif //TEMPO_DRPOLICY_HPP
