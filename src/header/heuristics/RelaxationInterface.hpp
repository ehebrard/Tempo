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
        ~AssumptionProxy() = default;

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
            s.assumption_level = s.numLiteral() - 1;
        }

    };

    /**
     * @brief Wrapper around AssumptionProxy that logs policy assumptions
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class LoggingAssumptionProxy {
        AssumptionProxy<T> proxy;
        std::vector<Literal<T>> assumptions{};
        std::filesystem::path logFile;
    public:
        /**
         * Ctor
         * @param logFile destination file
         * @param solver solver instance
         */
        LoggingAssumptionProxy(std::filesystem::path logFile, Solver<T> &solver) noexcept:
                proxy(solver), logFile(std::move(logFile)) {}

        /**
         * @copydoc AssumptionProxy::reset
         */
        void reset() {
            assumptions.clear();
            proxy.reset();
        }

        /**
         * @copydoc AssumptionProxy::makeAssumptions
         */
        template<concepts::typed_range<Literal<T>> L>
        bool makeAssumptions(L &&literals) {
            std::ranges::copy(std::forward<L>(literals), std::back_inserter(assumptions));
            return proxy.makeAssumptions(std::forward<L>(literals));
        }

        /**
         * @copydoc AssumptionProxy::tryMakeAssumption
         */
        bool tryMakeAssumption(Literal<T> lit) {
            assumptions.emplace_back(lit);
            return proxy.tryMakeAssumption(lit);
        }

        /**
         * @copydoc AssumptionProxy::getState
         * @note this triggers the serialization
         */
        [[nodiscard]] AssumptionState getState() const {
            serialization::serializeToFile(assumptions, logFile, std::ios_base::app);
            return proxy.getState();
        }
    };

    /**
     * @brief Assumption interface
     * @tparam I interface type
     * @tparam T timing type
     */
    template<typename I, typename T>
    concept AssumptionInterface = requires(I interface, Literal<T> l, std::vector<Literal<T>> lits) {
        interface.reset();
        { interface.makeAssumptions(lits) } -> std::convertible_to<bool>;
        { interface.tryMakeAssumption(l) } -> std::convertible_to<bool>;
        { interface.getState() } -> std::same_as<AssumptionState>;
    };

    /**
     * @brief Relaxation policy interface
     * @tparam P policy type
     * @tparam T timing type
     */
    template<typename P, typename T>
    concept RelaxationPolicy = requires(P policy, AssumptionProxy<T> interface, unsigned fails) {
        policy.relax(interface);
        policy.notifySuccess(fails);
        policy.notifyFailure();
    };
}

#endif //TEMPO_RELAXATIONINTERFACE_HPP
