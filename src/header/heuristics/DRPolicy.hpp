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
        instance.relax(proxy);
        { instance.newRegion() } -> std::convertible_to<bool>;
    };

    /**
     * @brief Repair policy interface
     */
    template<typename T>
    concept repair_policy = requires(T instance, unsigned numFail, AssumptionProxy<int> &proxy) {
        instance.notifySuccess(numFail);
        instance.notifyFailure(numFail);
        instance.notifyNewRegion();
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
        D destroy;
        R repair;
    public:

        /**
         * Ctor
         * @tparam ArgD destroy policy type
         * @tparam ArgR repair policy type
         * @param destroyer destroy policy
         * @param repairer repair policy
         */
        template<destroy_policy ArgD, repair_policy ArgR>
        RDPolicy(ArgD &&destroyer, ArgR &&repairer): destroy(std::forward<ArgD>(destroyer)),
                                                     repair(std::forward<ArgR>(repairer)) {}

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
            if (repair.exhausted()) {
                destroy.requestNewRegion();
            }

            destroy.relax(s);
            if (s.getState() == AssumptionState::Fail) {
                return;
            }

            if (destroy.newRegion()) {
                repair.notifyNewRegion();
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
     * @return constructed DR-policy
     */
    template<destroy_policy D, repair_policy R>
    auto make_RD_policy(D &&destroy, R &&repair) {
        return RDPolicy<D, R>(std::forward<D>(destroy), std::forward<R>(repair));
    }

    /**
     * @brief Random subset destroy policy.
     * @details @copybrief
     * Behaves like RandomSubset but only chooses a new region when asked to by the repair policy
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class RandomSubsetDestroy : protected RandomSubset<T> {
        bool exploreNewRegion = false;
    public:
        using RandomSubset<T>::notifySuccess;
        using RandomSubset<T>::notifyFailure;

        /**
         * Ctor
         * @param solver reference to solver
         * @param variables vector of boolean variables to choose from
         * @param ratio percentage of variables to destroy
         */
        RandomSubsetDestroy(Solver<T> &solver, std::vector<BooleanVar<T>> variables, double ratio) : RandomSubset<T>(
                solver, std::move(variables), ratio, 1.0) {}

        /**
         * requests the policy to choose a new region
         */
        void requestNewRegion() {
            this->fixed.clear();
        }

        /**
         * Destroy policy interface. Relaxes variables proportional to current destroy ration
         * @tparam AI assumption interface
         * @param proxy solver proxy
         */
        template<assumption_interface AI>
        void relax(AI &proxy) {
            if (this->fixed.empty()) {
                static_cast<RandomSubset<T>&>(*this).relax(proxy);
                exploreNewRegion = true;
            } else {
                proxy.makeAssumptions(this->fixed);
                exploreNewRegion = false;
            }
        }

        /**
         * Indicates whether a new region has been chosen
         * @return true if a new region was chosen in the last call to relax, false otherwise
         */
        [[nodiscard]] bool newRegion() const noexcept {
            return exploreNewRegion;
        }
    };
}

#endif //TEMPO_DRPOLICY_HPP
