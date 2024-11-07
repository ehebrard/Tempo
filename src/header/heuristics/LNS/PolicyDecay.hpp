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
                          bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit) noexcept;

        double fixRatio, decay, minFailRatio, maxFailRatio;
        bool decreaseOnSuccess;
        unsigned retryLimit;
        DecayMode decayMode;
    };

    std::ostream &operator<<(std::ostream &os, const PolicyDecayConfig &config);

    class PolicyDecay {
        PolicyDecayConfig config;
        unsigned failCount = 0;
        unsigned solverFailCount = 0;
        unsigned numLiterals;
        unsigned verbosity;
        double currentFixRatio;
    public:
        explicit PolicyDecay(const PolicyDecayConfig &config, unsigned numLiterals, unsigned verbosity) noexcept;

        [[nodiscard]] double decayFactor(double failRatio) const;

        [[nodiscard]] double calcFailRatio(unsigned fails) const noexcept;

        [[nodiscard]] unsigned getFailCount() const noexcept;

        void resetFailCount() noexcept;

        [[nodiscard]] double getFixRatio() const noexcept;

        void setFixRatio(double fixRatio) noexcept;

        void setNumLiterals(unsigned numLiterals) noexcept;

        void resetFixRatio() noexcept;

        void notifySuccess(unsigned numFails);

        void notifyFailure(unsigned numFails);

    };
}

#endif //TEMPO_POLICYDECAY_HPP
