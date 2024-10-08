/**
* @author Tim Luchterhand
* @date 04.10.24
* @brief
*/

#include "nn/GNNBackbonePredictor.hpp"

namespace tempo::nn {

    PolicyConfig::PolicyConfig(double relaxationRatio, double reactivity, double minCertainty, double minFailRatio,
                               double maxFailRatio, double confidenceUpdateFailRatio, bool carefulAssumptions,
                               bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit) noexcept
            : relaxationRatio(relaxationRatio), reactivity(reactivity), minCertainty(minCertainty),
              minFailRatio(minFailRatio), maxFailRatio(maxFailRatio),
              confidenceUpdateFailRatio(confidenceUpdateFailRatio), carefulAssumptions(carefulAssumptions),
              decreaseOnSuccess(decreaseOnSuccess), retryLimit(retryLimit), decayMode(decayMode) {}

    PolicyConfig::PolicyConfig() noexcept: relaxationRatio(1), reactivity(0.5), minCertainty(0.5),
                                           minFailRatio(-1), maxFailRatio(std::numeric_limits<double>::infinity()),
                                           confidenceUpdateFailRatio(std::numeric_limits<double>::infinity()),
                                           carefulAssumptions(false), decreaseOnSuccess(false), retryLimit(0),
                                           decayMode(DecayMode::Constant) {}

    std::ostream &operator<<(std::ostream &os, const PolicyConfig &config) {
        os << "-- GNN backbone predictor config:\n";
        os << "\t-- base relaxation ratio: " << config.relaxationRatio << "\n";
        os << "\t-- reactivity: " << config.reactivity << "\n";
        os << "\t-- GNN certainty threshold: " << config.minCertainty << "\n";
        os << "\t-- fail ratio interval: [" << config.minFailRatio << ", " << config.maxFailRatio << "]" << "\n";
        os << "\t-- GNN certainty update threshold: " << config.confidenceUpdateFailRatio << "\n";
        os << "\t-- careful assumptions on fail: " << (config.carefulAssumptions ? "true" : "false") << "\n";
        os << "\t-- decrease fix ration even on success: " << (config.decreaseOnSuccess ? "true" : "false") << "\n";
        os << "\t-- retry limit: " << config.retryLimit << "\n";
        os << "\t-- ratio decay mode: " << config.decayMode;
        return os;
    }
}