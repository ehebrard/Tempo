/**
 * @author Tim Luchterhand
 * @date 06.05.24
 * @brief Solution guided value selection
 */

#ifndef TEMPO_SOLUTIONGUIDED_HPP
#define TEMPO_SOLUTIONGUIDED_HPP

#include <vector>

#include "ReversibleObject.hpp"
#include "BaseBooleanHeuristic.hpp"
#include "heuristics/heuristic_interface.hpp"
#include "heuristics/TightestValue.hpp"
#include "heuristics/RandomBinaryValue.hpp"
#include "ReversibleObject.hpp"

namespace tempo::heuristics {

    namespace detail {
        template<typename Sched>
        concept solution_provider = requires(const Sched &s, var_t x) {
            { s.boolean.hasSolution() } -> std::convertible_to<bool>;
            { s.boolean.bestSolution() } -> std::same_as<const std::vector<bool> &>;
        };
    }

    /**
     * @brief Solution guided value selection heuristic.
     * @details @copybrief
     * Performs the first decent using some base heuristic. After that follows the most recent solution
     */
    template<class BaseHeuristic, concepts::scalar T>
    class SolutionGuided : public BaseBooleanHeuristic<SolutionGuided<BaseHeuristic, T>, T> {
        BaseHeuristic h;
        Reversible<size_t> var_ptr;
        Reversible<size_t> discrepancies;
        
        size_t max_discrepancies{2};
    public:
        template<class... Args>
        explicit SolutionGuided(Solver<T> &solver, Args &&...args) :
        BaseBooleanHeuristic<SolutionGuided, T>(solver.getOptions().polarity_epsilon), h(std::forward<Args>(args)...)
        , var_ptr(0, &(solver.getEnv()))
        , discrepancies(0, &(solver.getEnv()))
        {
            max_discrepancies = static_cast<size_t>(
                static_cast<double>(solver.boolean.size()) * solver.getOptions().sgd_ratio);
        }
        
        
        bool checkDiscrepancies(const Solver<T> &solver) const {
            const auto &b = solver.boolean;
            if (not b.hasSolution()) {
                return true;
            }
            
            const auto &sol = b.bestSolution();
            size_t num_discrepancies{0};
            auto n{static_cast<var_t>(solver.boolean.size())};
            for(var_t v{1}; v<n; ++v) {
                if(not solver.getBranch().has(v)) {
                    if(b.value(v) != sol[b.getLiteral(true, v)]) {
                        ++num_discrepancies;
                    }
                }
            }
            
            if(num_discrepancies != discrepancies) {
                std::cout << num_discrepancies << std::endl;
                for(var_t v{1}; v<n; ++v) {
                    std::cout << sol[b.getLiteral(true, v)] << " " ;
                    if(not solver.getBranch().has(v)) {
                        std::cout << b.value(v) ;
                        if(b.value(v) != sol[b.getLiteral(true, v)]) {
                            std::cout << " (" << v << ")" ;
                        }
                        std::cout << std::endl;
                    } else {
                        std::cout << "?\n";
                    }
                }
            }
            
            
            return num_discrepancies == discrepancies;
        }
        
        

        /**
         * heuristic interface
         * @tparam S class that provides previously encountered solutions
         */
        template<detail::solution_provider S>
        requires(heuristics::value_heuristic<BaseHeuristic, S> and boolean_info_provider<S>)
        [[nodiscard]] auto choose(var_t x, const S &solver) {
            const auto &b = solver.boolean;
            if (not b.hasSolution()) {
                return h.valueDecision({x, VariableType::Boolean}, solver);
            }

            auto& vars{solver.getBranch()};
            size_t cur_ptr{vars.start_idx()};
            
#ifdef DBG_DISCREPANCIES
            std::cout << std::endl << var_ptr << "/" << cur_ptr << " @" << solver.level() << std::endl;
#endif
            const auto &sol = b.bestSolution();
            
            auto vptr{static_cast<size_t>(var_ptr)};
            while(vptr < cur_ptr) {
                auto v{vars[vptr]};
                if(b.value(v) != sol[b.getLiteral(true, v)]) {
                    ++discrepancies;
#ifdef DBG_DISCREPANCIES
                    std::cout << " +" << v;
#endif
                }
#ifdef DBG_DISCREPANCIES
                else {
                    std::cout << " [" << v << "]";
                }
#endif
                ++vptr;
            }
            
//            std::cout << "\ndiscrepancies = " << discrepancies << std::endl;
            
            var_ptr = vptr; // the next decision is necessarily "correct"
            
            
#ifdef DBG_DISCREPANCIES
            if(not checkDiscrepancies(solver)) {
                std::cout << "bug!\n";
                exit(1);
            }
#endif
            
            assert(checkDiscrepancies(solver));
            
            auto posLit = b.getLiteral(true, x);
            auto negLit = b.getLiteral(false, x);
            
            
            if(sol.size() <= posLit or discrepancies > max_discrepancies) {
                return h.valueDecision({x, VariableType::Boolean}, solver);
            }
            
            assert(sol[posLit] != sol[negLit]);
            return sol[posLit] ? posLit : negLit;
        }
    };

    struct TSGFactory : MakeValueHeuristicFactory<TSGFactory> {
        TSGFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<SolutionGuided<TightestValue<T>, T>>(solver, solver);
        }
    };

    struct RSGFactory : MakeValueHeuristicFactory<RSGFactory> {
        RSGFactory();

        template<concepts::scalar T>
        [[nodiscard]] auto build_impl(Solver<T> &solver) const -> ValueHeuristic<T> {
            return std::make_unique<SolutionGuided<RandomBinaryValue<T>, T>>(solver);
        }
    };
}

#endif // TEMPO_SOLUTIONGUIDED_HPP
