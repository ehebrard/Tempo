/**
* @author Tim Luchterhand
* @date 07.11.24
* @brief
*/

#include <iostream>

#include "util/Options.hpp"
#include "heuristics/LNS/PolicyDecay.hpp"

namespace tempo::lns {
    PolicyDecayConfig::PolicyDecayConfig(double fixRatio, double baseDecay, double minFailRatio, double maxFailRatio,
                                         bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit,
                                         bool monotone)
        : fixRatio(fixRatio), decay(baseDecay), minFailRatio(minFailRatio), maxFailRatio(maxFailRatio),
          decreaseOnSuccess(decreaseOnSuccess), monotone(monotone), retryLimit(retryLimit), decayMode(decayMode) {
        if (fixRatio < 0 or fixRatio > 1) {
            throw std::invalid_argument("fixRatio must be between 0 and 1");
        }

        if (baseDecay < 0 or baseDecay > 1) {
            throw std::invalid_argument("baseDecay must be between 0 and 1");
        }
    }

    PolicyDecayConfig::PolicyDecayConfig() noexcept: fixRatio(1), decay(0.5), minFailRatio(-1),
                                                     maxFailRatio(std::numeric_limits<double>::infinity()),
                                                     decreaseOnSuccess(false), monotone(false), retryLimit(0),
                                                     decayMode(DecayMode::Constant) {}

    std::ostream &operator<<(std::ostream &os, const PolicyDecayConfig &config) {
        os << std::boolalpha;
        os << "-- policy decay config:\n";
        os << "\t-- ratio decay mode: " << config.decayMode << "\n";
        os << "\t-- monotone decent: " << config.monotone << "\n";
        os << "\t-- base fix ratio: " << config.fixRatio << "\n";
        os << "\t-- decay: " << config.decay << "\n";
        os << "\t-- fail ratio interval: [" << config.minFailRatio << ", " << config.maxFailRatio << "]" << "\n";
        os << "\t-- decrease fix ratio even on success: " << config.decreaseOnSuccess << "\n";
        os << "\t-- retry limit: " << config.retryLimit << "\n";
        return os;
    }

    PolicyDecay::PolicyDecay(const PolicyDecayConfig &config, unsigned numLiterals, unsigned verbosity) noexcept:
            config(config), numLiterals(numLiterals), verbosity(verbosity), currentFixRatio(config.fixRatio) {}

    double PolicyDecay::decayFactor(double failRatio) const {
        using enum lns::DecayMode;
        switch(config.decayMode) {
            case Constant:
                return config.decay;
            case Reciprog:
                return std::min(config.decay / failRatio * config.maxFailRatio, config.decay);
            case Exponential:
                return std::min(std::pow(config.decay, failRatio / config.maxFailRatio), config.decay);
            default:
                throw std::runtime_error("enum out of bounds");
        }
    }

    double PolicyDecay::calcFailRatio(unsigned int fails) const noexcept {
        return static_cast<double>(fails - solverFailCount) / numLiterals;
    }

    unsigned PolicyDecay::getFailCount() const noexcept {
        return failCount;
    }

    void PolicyDecay::resetFailCount() noexcept {
        failCount = 0;
    }

    double PolicyDecay::getFixRatio() const noexcept {
        return currentFixRatio;
    }

    void PolicyDecay::setFixRatio(double fixRatio) noexcept {
        currentFixRatio = fixRatio;
    }

    void PolicyDecay::resetFixRatio() noexcept {
        if (not config.monotone) {
            currentFixRatio = config.fixRatio;
        }
    }

    void PolicyDecay::notifyFailure(unsigned int numFails) {
        if (++failCount > config.retryLimit) {
            currentFixRatio *= decayFactor(calcFailRatio(numFails));
            if (verbosity >= Options::YACKING) {
                std::cout << std::setprecision(2) << "-- setting fix ratio = "
                          << currentFixRatio * 100 << "%" << std::endl;
            }
        } else {
            if (verbosity >= Options::YACKING) {
                std::cout << "-- backbone prediction failed, retrying more carefully\n";
            }
        }

        solverFailCount = numFails;
    }

    void PolicyDecay::notifySuccess(unsigned int numFails) {
        const auto failRatio = calcFailRatio(numFails);
        if (failRatio > config.maxFailRatio and config.decreaseOnSuccess) {
            currentFixRatio *= decayFactor(failRatio);
            if (verbosity >= Options::YACKING) {
                std::cout << "-- decreasing fix ratio to " << currentFixRatio * 100
                          << "% after too many solver fails" << std::endl;
            }
        } else if (failRatio < config.minFailRatio) {
            currentFixRatio = std::min(1.0, currentFixRatio / config.decay);
            if (verbosity >= Options::YACKING) {
                std::cout << "-- increasing fix ratio to " << currentFixRatio * 100
                          << "%" << std::endl;
            }
        }

        solverFailCount = numFails;
        failCount = 0;
    }

    void PolicyDecay::setNumLiterals(unsigned int numLiterals) noexcept {
        this->numLiterals = numLiterals;
    }
}