/**
* @author Tim Luchterhand
* @date 07.11.24
* @brief
*/

#ifndef TEMPO_POLICYDECAY_HPP
#define TEMPO_POLICYDECAY_HPP

#include <ostream>

#include "util/enum.hpp"
#include "heuristics/LNS/relaxation_interface.hpp"

namespace tempo::lns {

    PENUM(DecayMode, Constant, Reciprog, Exponential);

    struct PolicyDecayConfig {
        PolicyDecayConfig() noexcept;

        /**
         * Ctor
         * @param fixRatio percentage of literals to fix [0, 1]
         * @param baseDecay factor to apply to relaxation ratio on fail
         * @param minFailRatio lower bound solver failure rate at which to increase relaxation ratio
         * @param maxFailRatio upper bound solver failure rate at which to decrease relaxation ratio
         * @param decreaseOnSuccess whether to decrease fix rate even on success
         * @param retryLimit number of retries with same relaxation ration before decreasing relaxation ratio
         * @param decayMode type of decay to apply on fail or after too many solver fails
         */
        PolicyDecayConfig(double fixRatio, double baseDecay, double minFailRatio, double maxFailRatio,
                          bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit, bool monotone) noexcept;

        double fixRatio, decay, minFailRatio, maxFailRatio;
        bool decreaseOnSuccess, monotone;
        unsigned retryLimit;
        DecayMode decayMode;
    };

    std::ostream &operator<<(std::ostream &os, const PolicyDecayConfig &config);

    /**
     * @brief Class that decays a ratio value on different conditions using a specified decay scheme.
     */
    class PolicyDecay {
        PolicyDecayConfig config;
        unsigned failCount = 0;
        unsigned solverFailCount = 0;
        unsigned numLiterals;
        unsigned verbosity;
        double currentFixRatio;
    public:
        /**
         * Ctor
         * @param config decay configuration
         * @param numLiterals number of literals
         * @param verbosity verbosity for logging
         */
        explicit PolicyDecay(const PolicyDecayConfig &config, unsigned numLiterals, unsigned verbosity) noexcept;

        /**
         * Calculate the decay factor given a fail ratio according to the decay policy
         * @param failRatio fail ratio
         * @return the decay factor according to the decay policy
         */
        [[nodiscard]] double decayFactor(double failRatio) const;

        /**
         * Calculate the fail ratio (number of fails since the last iterations with respect to the number of literals)
         * @param fails current number of solver fails
         * @return fail ratio
         */
        [[nodiscard]] double calcFailRatio(unsigned fails) const noexcept;

        /**
         * Get the current LNS iterations fail count
         * @return
         */
        [[nodiscard]] unsigned getFailCount() const noexcept;

        /**
         * Manually reset the failed LNS iterations count (not the solver fail count)
         */
        void resetFailCount() noexcept;

        /**
         * Get the current fix ratio
         * @return current fix ratio
         */
        [[nodiscard]] double getFixRatio() const noexcept;

        /**
         * Manually set the fix ratio
         * @param fixRatio new value
         */
        void setFixRatio(double fixRatio) noexcept;

        /**
         * Manually set the number of literals
         * @param numLiterals
         */
        void setNumLiterals(unsigned numLiterals) noexcept;

        /**
         * Reset the fix ratio to the value specified in the config
         */
        void resetFixRatio() noexcept;

        /**
         * Call this method after a successful LNS iteration to update the fix ratio according to the decay policy
         * @param numFails number of solver fails thus far
         */
        void notifySuccess(unsigned numFails);

        /**
         * Call this method after a failed LNS iteration to update the fix ratio according to the decay policy
         * @param numFails number of solver fails thus far
         */
        void notifyFailure(unsigned numFails);

    };
}

#endif //TEMPO_POLICYDECAY_HPP
