/**
* @author Tim Luchterhand
* @date 16.12.24
* @file PerfectRelaxationOracle.hpp
* @brief Contains the perfect relaxation heuristic
*/

#ifndef PERFECTRELAXATIONORACLE_HPP
#define PERFECTRELAXATIONORACLE_HPP

#include "util/traits.hpp"
#include "util/random.hpp"
#include "Solution.hpp"
#include "heuristics/LNS/relaxation_interface.hpp"


namespace tempo::lns {
    /**
     * @brief Perfect relaxation heuristic. Follows a target solution and relaxes only edges that are inconsistent
     * with said solution
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class PerfectRelaxationOracle {
        Solution<T> solution;
        std::vector<BooleanVar<T>> variables;
        NumericVar<T> schedule;
        double fixRatio;
        double epsilon;
    public:
        /**
         * Ctor
         * @param solution the solution to follow
         * @param variables binary variables in the problem
         * @param scheduleDuration schedule duration variable
         * @param fixRatio percentage of variables to fix
         * @param epsilon error probability. 0 means that the oracle will never
         * fix edges inconsistent with the target solution
         */
        PerfectRelaxationOracle(Solution<T> solution, NumericVar<T> scheduleDuration,
                                std::vector<BooleanVar<T>> variables, double fixRatio,
                                double epsilon) : solution(std::move(solution)),
                                                  variables(std::move(variables)), schedule(scheduleDuration),
                                                  fixRatio(fixRatio), epsilon(epsilon) {
            if (epsilon < 0 or epsilon > 1) {
                throw std::invalid_argument("epsilon must be between 0 and 1");
            }

            if (fixRatio < 0 or fixRatio > 1) {
                throw std::invalid_argument("fixRatio must be between 0 and 1");
            }
        }

        PerfectRelaxationOracle(const PerfectRelaxationOracle&) = default;
        PerfectRelaxationOracle(PerfectRelaxationOracle&&) = default;
        PerfectRelaxationOracle &operator=(const PerfectRelaxationOracle&) = default;
        PerfectRelaxationOracle &operator=(PerfectRelaxationOracle&&) = default;

        /**
         * Relaxation policy interface
         * @tparam AI assumption proxy type
         * @param proxy
         */
        template<assumption_interface AI>
        void relax(AI &proxy) {
            using namespace std::views;
            if (proxy.getSolver().numeric.hasSolution() and proxy.getSolver().numeric.lower(schedule) == solution.
                numeric.lower(schedule)) {
                if (proxy.getSolver().getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- oracle: optimum reached" << std::endl;
                }
                return;
            }
            const auto numFix = static_cast<std::size_t>(variables.size() * fixRatio);
            std::ranges::shuffle(variables, RNG{});
            auto selection = variables | filter([this, &proxy](auto bv) {
                return solution.boolean.value(bv) == proxy.getSolver().boolean.value(bv) or
                       random_event_occurred(epsilon);
            }) | take(numFix) | transform([&proxy](auto bv) {
                return bv == proxy.getSolver().boolean.value(bv);
            });

            proxy.makeAssumptions(selection);
        }

        void notifySuccess(unsigned) noexcept {}

        void notifyFailure(unsigned) noexcept {}
    };

}

#endif //PERFECTRELAXATIONORACLE_HPP
