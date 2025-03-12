/**
* @author Tim Luchterhand
* @date 11.03.25
* @file GNNValueHeuristic.cpp
* @brief
*/

#include "nn/GNNValueHeuristics.hpp"


namespace tempo::nn {
    std::ostream &operator<<(std::ostream &os, const DispatcherConfig &config) {
        auto precision = os.precision();
        os << std::setprecision(2);
        os << "-- Dispatcher configuration:" << std::endl;
        os << "\t-- Ratio increment on fail: " << config.failIncrement * 100 << "% trigger threshold" << std::endl;
        os << "\t-- Ratio decrement on success: " << config.successDecrement * 100 << "% of trigger threshold" << std::endl;
        os << "\t-- Ratio increment on restart: " << config.restartIncrement * 100 << "% of trigger threshold" << std::endl;
        os << "\t-- Ratio decrement on solution: " << config.solutionDecrement * 100 << "% of trigger threshold" <<
                std::endl;
        os << "\t-- Max ratio: " << config.maxFillRate * 100 << "% of trigger threshold" << std::endl;
        os << "\t-- Heat increment: " << config.heatIncrement * 100 << "% of max heat" << std::endl;
        os << "\t-- Heat lower threshold: " << config.heatLowerThreshold* 100 << "% of max heat" << std::endl;
        os << "\t-- Heat decay: " << config.heatDecay;
        os << std::setprecision(static_cast<int>(precision));
        return os;
    }

    Dispatcher::Dispatcher(const DispatcherConfig &config) noexcept
        : failIncrement(static_cast<int>(TriggerThreshold * config.failIncrement)),
          successDecrement(static_cast<int>(TriggerThreshold * config.successDecrement)),
          restartIncrement(static_cast<int>(TriggerThreshold * config.restartIncrement)),
          solutionDecrement(static_cast<int>(TriggerThreshold * config.solutionDecrement)),
          maxFillRate(static_cast<int>(TriggerThreshold * config.maxFillRate)),
          heatIncrement(config.heatIncrement), heatDecay(config.heatDecay),
          heatLowerThreshold(config.heatLowerThreshold) {}

    void Dispatcher::onFail() noexcept {
        fillingRate = std::clamp(fillingRate + failIncrement, -maxFillRate, maxFillRate);
    }

    void Dispatcher::onSuccess() noexcept {
        fillingRate = std::clamp(fillingRate - successDecrement, -maxFillRate, maxFillRate);
    }

    void Dispatcher::onRestart() noexcept {
        fillingRate = std::clamp(fillingRate + restartIncrement, -maxFillRate, maxFillRate);
    }

    void Dispatcher::onSolution() noexcept {
        fillingRate = std::clamp(fillingRate - solutionDecrement, -maxFillRate, maxFillRate);
    }

    void Dispatcher::onInference() noexcept {
        waterLevel = 0;
        temperature += heatIncrement;
        if (temperature >= MaxTemperature) {
            isOverheated = true;
        }
    }

    void Dispatcher::step() noexcept {
        if (fillingRate > 0 or std::abs(fillingRate) <= waterLevel) {
            waterLevel += fillingRate;
        }

        temperature *= heatDecay;
        if (temperature < heatLowerThreshold) {
            isOverheated = false;
        }
    }

    bool Dispatcher::inferenceAllowed() const noexcept {
        return not isOverheated and waterLevel >= TriggerThreshold;
    }

    bool Dispatcher::overheated() const noexcept {
        return isOverheated;
    }

    std::ostream & operator<<(std::ostream &os, const Dispatcher &dispatcher) {
        auto precision = os.precision();
        os << std::setprecision(2);
        os << "-- GNN dispatcher info:" << std::endl;
        os << "\t-- fill level: " << 100 * static_cast<double>(dispatcher.waterLevel) / Dispatcher::TriggerThreshold <<
                "%" << std::endl;
        os << "\t-- fill rate: " << 100 * static_cast<double>(dispatcher.fillingRate) / Dispatcher::TriggerThreshold <<
                "%/step" << std::endl;
        os << "\t-- heat: " << 100 * static_cast<double>(dispatcher.temperature) / Dispatcher::MaxTemperature << "%";
        os << std::setprecision(static_cast<int>(precision));
        return os;
    }
}
