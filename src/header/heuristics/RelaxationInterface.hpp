/**
* @author Tim Luchterhand
* @date 19.09.24
* @brief LNS relaxation policy interface.
*/

#ifndef TEMPO_RELAXATIONINTERFACE_HPP
#define TEMPO_RELAXATIONINTERFACE_HPP

#include <concepts>
#include <filesystem>
#include <vector>

#include "util/traits.hpp"
#include "util/serialization.hpp"
#include "Literal.hpp"
#include "Failure.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::heuristics {

    /**
     * @brief State of the AssumptionInterface used by the solver to determine followup actions
     */
    enum class AssumptionState {
        Success, ///< relaxation succeeded
        Fail, ///< relaxation failed and relaxation policy cannot handle the error
        Empty ///< relaxation policy does not want to make assumptions
    };

    /**
     * @brief Assumption interface
     * @tparam I interface type
     */
    template<typename I>
    concept assumption_interface = requires(I interface, Literal<int> l, std::vector<Literal<int>> lits) {
        interface.reset();
        { interface.makeAssumptions(lits) } -> std::convertible_to<bool>;
        { interface.tryMakeAssumption(l) } -> std::convertible_to<bool>;
        { interface.getState() } -> std::same_as<AssumptionState>;
        { interface.getSolver() } -> concepts::same_template_clvref<Solver>;
    };

    /**
     * @brief Wrapper around solver that only exposes methods necessary for making assumptions.
     * @detail @copybrief
     * Ensures propagation and correct assumption level at all times
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class AssumptionProxy {
        Solver<T> &s;
        int baseLevel;
        AssumptionState state = AssumptionState::Empty;

    public:

        AssumptionProxy(const AssumptionProxy &) = delete;
        AssumptionProxy(AssumptionProxy &&) = delete;
        AssumptionProxy &operator=(const AssumptionProxy &) = delete;
        AssumptionProxy &operator=(AssumptionProxy &&) = delete;
        ~AssumptionProxy() { s.assumption_stamp = s.ground_stamp; };

        /**
         * Ctor
         * @param solver
         */
        AssumptionProxy(Solver<T> &solver): s(solver) {
            s.initializeSearch();
            baseLevel = s.saveState();
        }

        /**
         * Undoes all assumptions made to this point.
         * @note call this function to get out of inconsistent state
         */
        void reset() {
            stateTransition(AssumptionState::Empty);
            s.restoreState(baseLevel);
        }

        /**
         * Sets all literals in the given range
         * @tparam L literal range type
         * @param literals range of literals
         * @return true if all literals were set with success, false otherwise
         * @throws std::runtime_error if in inconsistent state
         * @note in order to retry and make other assumptions, reset() needs to be called
         */
        template<concepts::typed_range<Literal<T>> L>
        bool makeAssumptions(L &&literals) {
            using enum AssumptionState;
            if (state == Fail) {
                throw std::runtime_error("cannot make assumptions in inconsistent state. Did you forget to try reset?");
            }

            if (std::ranges::empty(std::forward<L>(literals))) {
                return true;
            }

            try {
                for (auto &&lit : std::forward<L>(literals)) {
                    s.set(std::forward<decltype(lit)>(lit));
                }

                s.propagate();
            } catch (const Failure<T> &) {
                stateTransition(Fail);
                return false;
            }

            stateTransition(Success);
            setAssumptionLevel();
            return true;
        }

        /**
         * Tries to set a single literal
         * @param lit literal to assume
         * @return true if set was successful, false if inconsistent
         * @throws std::runtime_error if in inconsistent state
         */
        bool tryMakeAssumption(Literal<T> lit) {
            using enum AssumptionState;
            if (state == Fail) {
                throw std::runtime_error("cannot make assumptions in inconsistent state. Did you forget to try reset?");
            }

            if ((not lit.isBoolean() or s.boolean.falsified(lit)) and (lit.isBoolean() or s.numeric.falsified(lit))) {
                return false;
            }

            auto cp = s.saveState();
            try {
                s.post(lit);
            } catch(const Failure<T> &) {
                s.restoreState(cp);
                return false;
            }

            stateTransition(AssumptionState::Success);
            setAssumptionLevel();
            return true;
        }

        /**
         * Gets the current state of the assumptions
         * @return
         */
        [[nodiscard]] AssumptionState getState() const noexcept {
            return state;
        }

        /**
         * const reference to the solver
         * @return const reference to the solver
         */
        [[nodiscard]] const auto &getSolver() const noexcept {
            return s;
        }

    private:

        void stateTransition(AssumptionState as) noexcept {
            using enum AssumptionState;
            switch (state) {
                case Empty:
                    state = as;
                    break;
                case Fail:
                    state = as == Empty ? as : state;
                    break;
                case Success:
                    state = as;
            }
        }

        void setAssumptionLevel() noexcept {
          s.assumption_stamp = s.numLiteral();
        }

    };

    /**
     * @brief Assumption interface wrapper that collects the assumptions made
     * @tparam T timing type
     * @tparam AI assumption interface type
     */
    template<concepts::scalar T, assumption_interface AI>
    class AssumptionCollector {
        std::vector<Literal<T>> assumptions;
        AI &proxy;
    public:
        AssumptionCollector(const AssumptionCollector &) = delete;
        AssumptionCollector(AssumptionCollector &&) = delete;
        AssumptionCollector &operator=(const AssumptionCollector &) = delete;
        AssumptionCollector &operator=(AssumptionCollector &&) = delete;

        /**
         * Ctor
         * @param proxy assumption interface to wrap
         */
        AssumptionCollector(AI &proxy) noexcept: proxy(proxy) {}

        /**
         * Calls reset on base assumption interface. Also clears the assumption cache
         */
        void reset() {
            proxy.reset();
            assumptions.clear();
        }

        /**
         * Calls tryMakeAssumption on base assumption interface
         * @param lit literal to set
         * @return true if set was successful, false otherwise
         */
        bool tryMakeAssumption(Literal<T> lit) {
            bool res = proxy.tryMakeAssumption(lit);
            if (res) {
                assumptions.emplace_back(lit);
            }

            return res;
        }

        /**
         * Calls makeAssumptions on base assumption interface
         * @tparam L literal range type
         * @param literals range of literals to set
         * @return true if set was successful, false otherwise
         */
        template<concepts::typed_range<Literal<T>> L>
        bool makeAssumptions(L &&literals) {
            bool res = proxy.makeAssumptions(literals);
            if (res) {
                std::ranges::copy(std::forward<L>(literals), std::back_inserter(assumptions));
            }

            return res;
        }

        /**
         * Gets the current state of the assumptions
         * @return
         */
        [[nodiscard]] AssumptionState getState() const noexcept {
            return proxy.getState();
        }

        /**
         * Returns the assumptions made since the last call to rest
         * @return
         */
        auto getAssumptions() noexcept -> std::vector<Literal<T>> & {
            return assumptions;
        }

        /**
         * const reference to the solver
         * @return const reference to the solver
         */
        [[nodiscard]] decltype(auto) getSolver() const noexcept {
            return proxy.getSolver();
        }
    };

    enum class PolicyAction {
        Reset,
        TrySet,
        Set
    };

    /**
     * @brief Wrapper around AssumptionProxy that logs policy assumptions.
     * @tparam T timing type
     * @note This class is slower than the base AssumptionProxy because it flushes all decisions directly to a log file.
     * It should only be used for debugging. Also, due to the way LNS is implemented, this class gets created and
     * destroyed multiple times and is therefore unable to create valid json. You need to at least insert one
     * opening '[' ath the beginning and one closing ']' at the end of the log file. If the program crashed, then most
     * probably there are further missing ']' at the end of the file. But the rest should™ be valid json.
     */
    template<concepts::scalar T>
    class LoggingAssumptionProxy {
    public:
        using PolicyTrace = std::vector<std::pair<PolicyAction, std::vector<Literal<T>>>>;
    private:
        AssumptionProxy<T> proxy;
        std::ofstream logFile;
        bool firstAction = true;

        template<serialization::serializable S>
        void flushToLog(const S &object) {
            if (not firstAction) {
                logFile << ",";
            }

            firstAction = false;
            nlohmann::json j = object;
            logFile << j.dump(__JSON_INDENT__) << std::flush;
        }

    public:
        LoggingAssumptionProxy(const LoggingAssumptionProxy &) = default;
        LoggingAssumptionProxy(LoggingAssumptionProxy &&) = default;
        LoggingAssumptionProxy &operator=(const LoggingAssumptionProxy &) = default;
        LoggingAssumptionProxy &operator=(LoggingAssumptionProxy &&) = default;

        ~LoggingAssumptionProxy() {
            logFile << "],";
        }

        /**
         * Ctor
         * @param logFile destination file
         * @param solver solver instance
         */
        LoggingAssumptionProxy(const std::filesystem::path &logFile, Solver<T> &solver):
                proxy(solver), logFile(logFile, std::ios_base::app) {
            if (not this->logFile.is_open()) {
                throw std::runtime_error("could not open log file for policy logging");
            }

            this->logFile << "[";
        }

        /**
         * @copydoc AssumptionProxy::reset
         */
        void reset() {
            auto data = std::make_pair(PolicyAction::Reset, std::vector<Literal<T>>{});
            flushToLog(data);
            proxy.reset();
        }

        /**
         * @copydoc AssumptionProxy::makeAssumptions
         */
        template<concepts::typed_range<Literal<T>> L>
        bool makeAssumptions(L &&literals) {
            std::vector<Literal<T>> lits;
            lits.reserve(std::ranges::size(literals));
            std::ranges::copy(std::forward<L>(literals), std::back_inserter(lits));
            auto data = std::make_pair(PolicyAction::Set, std::move(lits));
            flushToLog(data);
            return proxy.makeAssumptions(std::forward<L>(literals));
        }

        /**
         * @copydoc AssumptionProxy::tryMakeAssumption
         */
        bool tryMakeAssumption(Literal<T> lit) {
            auto data = std::make_pair(PolicyAction::TrySet, std::vector{lit});
            flushToLog(data);
            return proxy.tryMakeAssumption(lit);
        }

        /**
         * @copydoc AssumptionProxy::getState
         * @note this triggers the serialization
         */
        [[nodiscard]] AssumptionState getState() const {
            return proxy.getState();
        }

        /**
         * const reference to the solver
         * @return const reference to the solver
         */
        [[nodiscard]] decltype(auto) getSolver() const noexcept {
            return proxy.getSolver();
        }
    };

    /**
     * @brief Relaxation policy interface
     * @tparam P policy type
     * @tparam T timing type
     */
    template<typename P>
    concept relaxation_policy = requires(P policy, AssumptionProxy<int> interface, unsigned fails) {
        policy.relax(interface);
        policy.notifySuccess(fails);
        policy.notifyFailure(fails);
    };
}

#endif //TEMPO_RELAXATIONINTERFACE_HPP
