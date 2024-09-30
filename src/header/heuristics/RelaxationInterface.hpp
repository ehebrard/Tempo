/**
* @author Tim Luchterhand
* @date 19.09.24
* @brief LNS relaxation policy interface.
*/

#ifndef TEMPO_RELAXATIONINTERFACE_HPP
#define TEMPO_RELAXATIONINTERFACE_HPP

#include "util/traits.hpp"
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
    template<typename T>
    class AssumptionInterface {
        Solver<T> &s;
        int baseLevel;
        AssumptionState state = AssumptionState::Empty;

    public:
        AssumptionInterface(const AssumptionInterface &) = delete;
        AssumptionInterface(AssumptionInterface &&) = delete;
        AssumptionInterface &operator=(const AssumptionInterface &) = delete;
        AssumptionInterface &operator=(AssumptionInterface &&) = delete;
        ~AssumptionInterface() = default;

        /**
         * Ctor
         * @param solver
         */
        AssumptionInterface(Solver<T> &solver): s(solver) {
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

            for (auto &&lit : std::forward<L>(literals)) {
                s.set(std::forward<decltype(lit)>(lit));
            }

            try {
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
     * @brief Relaxation policy interface
     * @tparam P policy type
     * @tparam T timing type
     */
    template<typename P, typename T>
    concept RelaxationPolicy = requires(P policy, AssumptionInterface<T> interface) {
        policy.relax(interface);
        policy.notifySuccess();
        policy.notifyFailure();
    };
}

#endif //TEMPO_RELAXATIONINTERFACE_HPP
