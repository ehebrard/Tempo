/**
* @author Tim Luchterhand
* @date 04.10.24
* @brief
*/

#include "nn/GNNBackbonePredictor.hpp"

namespace tempo::nn {

    PolicyConfig::PolicyConfig(double fixRatio, double reactivity, double minCertainty, double minFailRatio,
                               double maxFailRatio, double exhaustionThreshold, lns::AssumptionMode assumptionMode,
                               bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit) noexcept
            : fixRatio(fixRatio), decay(reactivity), minCertainty(minCertainty),
              minFailRatio(minFailRatio), maxFailRatio(maxFailRatio), exhaustionThreshold(exhaustionThreshold),
              decreaseOnSuccess(decreaseOnSuccess), retryLimit(retryLimit), decayMode(decayMode),
              assumptionMode(assumptionMode) {}

    PolicyConfig::PolicyConfig() noexcept: fixRatio(1), decay(0.5), minCertainty(0.5),
                                           minFailRatio(-1), maxFailRatio(std::numeric_limits<double>::infinity()),
                                           exhaustionThreshold(0.01), decreaseOnSuccess(false), retryLimit(0),
                                           decayMode(DecayMode::Constant), assumptionMode(lns::AssumptionMode::BestN) {}

    std::ostream &operator<<(std::ostream &os, const PolicyConfig &config) {
        os << "-- GNN backbone predictor config:\n";
        os << "\t-- base fix ratio: " << config.fixRatio << "\n";
        os << "\t-- decay: " << config.decay << "\n";
        os << "\t-- GNN certainty threshold: " << config.minCertainty << "\n";
        os << "\t-- fail ratio interval: [" << config.minFailRatio << ", " << config.maxFailRatio << "]" << "\n";
        os << "\t-- exhaustion threshold: " << config.exhaustionThreshold << "\n";
        os << "\t-- assumption mode: " << config.assumptionMode << "\n";
        os << "\t-- decrease fix ratio even on success: " << (config.decreaseOnSuccess ? "true" : "false") << "\n";
        os << "\t-- retry limit: " << config.retryLimit << "\n";
        os << "\t-- ratio decay mode: " << config.decayMode;
        return os;
    }
}