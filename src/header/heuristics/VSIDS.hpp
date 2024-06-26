//
// Created by tim on 15.11.22.
//

#ifndef TEMPO_VSIDS_HPP
#define TEMPO_VSIDS_HPP

#include "Global.hpp"
#include "RankingHeuristic.hpp"
#include "heuristics/impl/DecayingEventActivityMap.hpp"
#include "util/SubscribableEvent.hpp"

namespace tempo::heuristics {
    /**
     * @brief VSIDS (Variable State Independent Decaying Sum) heuristic
     * @details @copybrief
     * Selects variables which occur in the largest number of clauses (proportionally to highest activity). Activity for
     * a literal l is calculated the following way:
     * \f$ A_l = \sum_{i=1}^n \gamma^{n - i} w_i \f$ where \f$ \gamma \f$ is a constant decay factor and \f$ n \f$ is
     * the the number of calls to the function VSIDS::updateActivity
     * The final score of a literal is inversely proportional to its activity (@see getCost)
     */
    template<typename T>
    class VSIDS : public RankingHeuristic<VSIDS<T>> {
    public:


        /**
         * CTor. Infers parameters from given scheduler. Subscribes to events of interest
         * @param solver target solver
         */
        explicit VSIDS(Solver<T> &solver) :
                activity(solver, solver.getOptions().vsids_decay),
                handlerToken(solver.ClauseAdded.subscribe_handled(
                        [this, &solver](const auto &arg) { this->activity.update(arg, solver); })) {}

        // Copy and move are disabled. Otherwise, a call to the subscribed event handler will cause undefined behavior
        VSIDS(const VSIDS &) = delete;

        VSIDS(VSIDS &&) = delete;

        VSIDS &operator=(const VSIDS &) = delete;

        VSIDS &operator=(VSIDS &&) = delete;

        ~VSIDS() = default;

        [[nodiscard]] double getCost(const var_t x, const Solver<T> &solver) const {

            //@TODO: there shoud be a normalization thingy and Boolean variables without semantic should get the highest value
            double dom{1};

            if (solver.boolean.hasSemantic(x)) {
                auto p{solver.boolean.getLiteral(true, x)};
                auto n{solver.boolean.getLiteral(false, x)};

                auto prec_a{solver.boolean.getEdge(p)};
                auto prec_b{solver.boolean.getEdge(n)};

                auto gap_a = solver.numeric.upper(prec_a.from) - solver.numeric.lower(prec_a.to);
                auto gap_b = solver.numeric.upper(prec_b.from) - solver.numeric.lower(prec_b.to);

                dom = static_cast<double>(std::max(gap_a, gap_b));
            }
            auto act{activity.get(x, solver)};

            //            std::cout << x << ": " << dom << "/" << act << "=" <<
            //            (dom/act) << std::endl;

            return dom / act;
        }

        /**
         * @param solver
         * @todo currently only selects boolean variables
         */
        auto nextVariable(const Solver<T> &solver) const -> VariableSelection {
            return {this->bestVariable(solver.getBranch(), solver), VariableType::Boolean};
        }

    private:
        impl::DecayingEventActivityMap<T> activity;
        SubscriberHandle handlerToken;
    };
}

#endif //TEMPO_VSIDS_HPP
