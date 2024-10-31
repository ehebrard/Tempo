/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief GNN based relaxation policy based on edge predictor
*/

#ifndef TEMPO_GNNBACKBONEPREDICTOR_HPP
#define TEMPO_GNNBACKBONEPREDICTOR_HPP

#include <filesystem>
#include <vector>
#include <ranges>
#include <iostream>
#include <limits>
#include <variant>
#include <Iterators.hpp>

#include "heuristics/RelaxationInterface.hpp"
#include "util/traits.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "Literal.hpp"
#include "GNNPrecedencePredictor.hpp"
#include "util/Profiler.hpp"
#include "util/enum.hpp"
#include "util/Options.hpp"
#include "util/factory_pattern.hpp"
#include "Solver.hpp"
#include "Constant.hpp"

namespace tempo::nn {
    namespace fs = std::filesystem;

    PENUM(DecayMode, Constant, Reciprog, Exponential);
    PENUM(AssumptionMode, SingleShot, GreedySkip, GreedyInverse, Optimal)

    struct PolicyConfig {
        PolicyConfig() noexcept;

        /**
         * Ctor
         * @param fixRatio percentage of literals to fix [0, 1]
         * @param reactivity factor to apply to relaxation ratio on fail
         * @param minCertainty minimum GNN certainty [0, 1]
         * @param minFailRatio lower bound solver failure rate at which to increase relaxation ratio
         * @param maxFailRatio upper bound solver failure rate at which to decrease relaxation ratio
         * @param exhaustionThreshold lower bound on ratio of literals at which a new region should be explored
         * @param assumptionMode how to make assumptions
         * @param decreaseOnSuccess whether to decrease fix rate even on success
         * @param retryLimit number of retries with same relaxation ration before decreasing relaxation ratio
         * @param decayMode type of decay to apply on fail or after too many solver fails
         */
        PolicyConfig(double fixRatio, double reactivity, double minCertainty, double minFailRatio,
                     double maxFailRatio, double exhaustionThreshold, AssumptionMode assumptionMode,
                     bool decreaseOnSuccess, DecayMode decayMode, unsigned int retryLimit) noexcept;

        double fixRatio, decay, minCertainty, minFailRatio, maxFailRatio, exhaustionThreshold;
        bool decreaseOnSuccess;
        unsigned retryLimit;
        DecayMode decayMode;
        AssumptionMode assumptionMode;
    };

    std::ostream &operator<<(std::ostream &os, const PolicyConfig &config);

    /**
     * @brief Fix policy that simply takes the best N edges predicted by the GNN
     */
    struct BestN {
        constexpr void reset() noexcept {};

        template<concepts::scalar T, heuristics::assumption_interface AI>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned,
                           const std::vector<std::pair<Literal<T>, double>> & gnnOutput) {
            using namespace std::views;
            proxy.makeAssumptions(gnnOutput | take(numLiterals) | elements<0>);
            return std::min(numLiterals, gnnOutput.size());
        }
    };


    /**
     * @brief Fix policy that tries to fix N edges in a greedy fashion
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class GreedyFix {
        std::vector<Literal<T>> cache;
        bool inverse;
    public:
        /**
         * Ctor
         * @param inverse whether to set the inverse literal on fail
         */
        explicit constexpr GreedyFix(bool inverse) noexcept: inverse(inverse) {}

        void reset() noexcept {
            cache.clear();
        }

        template<heuristics::assumption_interface AI>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned numFails,
                    const std::vector<std::pair<Literal<T>, double>> & gnnOutput) {
            if (numFails == 0 and not cache.empty()) {
                proxy.makeAssumptions(cache | std::views::take(numLiterals));
                return std::min(numLiterals, cache.size());
            }

            if (proxy.getSolver().getOptions().verbosity >= Options::SOLVERINFO) {
                std::cout << "-- greedy selection of literals\n";
            }

            std::size_t litCount = 0;
            std::vector<Literal<T>> newAssumptions;
            newAssumptions.reserve(numLiterals);
            for (auto lit: gnnOutput | std::views::elements<0>) {
                if (litCount == numLiterals) {
                    break;
                }

                bool success = proxy.tryMakeAssumption(lit);
                if (not success and inverse) {
                    lit = ~lit;
                    success = proxy.tryMakeAssumption(lit);
                    if (not success) {
                        proxy.fail();
                        return 0;
                    }
                }

                litCount += success;
                if (success) {
                    newAssumptions.emplace_back(lit);
                }
            }

            std::swap(newAssumptions, cache);
            return litCount;
        }
    };

    namespace detail {
        template<concepts::scalar T, typename Lookup>
        class CachedLookup {
            std::vector<BooleanVar<T>> cache;
            Lookup lookup;
            Solver<T> &solver;
            const std::vector<BooleanVar<T>> &vars;
            static constexpr auto NoVar = BooleanVar<T>(Constant::NoVar, Constant::NoSemantic);
        public:
            template<typename L>
            CachedLookup(Solver<T> &solver, const std::vector<BooleanVar<T>> &vars,
                         std::size_t numLits,L &&lookup): lookup(std::forward<L>(lookup)), solver(solver), vars(vars) {
                assert(vars.size() == this->lookup.size());
                cache.resize(numLits, NoVar);
            }

            template<concepts::scalar U>
            auto operator()(Literal<U> lit) -> Literal<T> {
                auto &var = cache[lit.variable()];
                if (var != NoVar) {
                    return var == lit.sign();
                }

                auto res = std::ranges::find(lookup, lit);
                if (res != std::ranges::end(lookup)) {
                    var = vars[std::ranges::distance(std::ranges::begin(lookup), res)];
                } else {
                    //@TODO are some variables constants in the sub problem?
                    var = solver.newBoolean();
                    solver.addToSearch(var);
                }

                return var == lit.sign();
            }
        };

        template<concepts::scalar T, typename Lookup>
        auto makeLookup(Solver<T> &solver, const std::vector<BooleanVar<T>> &vars,
                        std::size_t numLits, Lookup &&lookup) {
            return CachedLookup<T, Lookup>(solver, vars, numLits, std::forward<Lookup>(lookup));
        }
    }

    /**
     * @brief Fix policy that tries to fix N edges by maximizing the sum of their confidence values while respecting
     * the learned clauses
     * @todo integrate problem constraints
     * @tparam T timing type
     */
    template<concepts::scalar T>
    class OptimalFix {
        std::vector<Literal<T>> cache;
        unsigned timeLimit;
        unsigned failLimit;
        bool hardTimeLimit;
        static constexpr auto FixPointPrec = 1000; //@TODO use Solver<float>
    public:

        /**
         * Ctor
         * @param timeLimit time limit in ms for optimization problem
         * @param failLimit fail limit for optimization problem
         * @param hardTimeLimit whether to always cut off after the time limit even when no solution has been found
         */
        OptimalFix(unsigned timeLimit, unsigned failLimit, bool hardTimeLimit) noexcept: timeLimit(timeLimit),
                                                                                         failLimit(failLimit),
                                                                                         hardTimeLimit(hardTimeLimit) {}

        void reset() noexcept {
            cache.clear();
        }

        template<heuristics::assumption_interface AI>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned numFails,
                    const std::vector<std::pair<Literal<T>, double>> & gnnOutput) {
            using namespace std::views;
            using ST = int;
            if (numFails == 0 and not cache.empty()) {
                proxy.makeAssumptions(cache | take(numLiterals));
                return std::min(numLiterals, cache.size());
            }

            cache.clear();
            Options opt;
            opt.search_limit = failLimit;
            opt.verbosity = Options::SILENT;
            if (proxy.getSolver().getOptions().verbosity >= Options::SOLVERINFO) {
                std::cout << "-- solving selection problem\n";
                opt.verbosity = Options::NORMAL;
            }

            tempo::util::StopWatch stopWatch;
            bool solution = false;
            Solver<ST> solver(opt);
            solver.SolutionFound.subscribe_unhandled([&solution, &solver, &stopWatch, this](auto &&) {
                solution = true;
                if (stopWatch.elapsed<std::chrono::milliseconds>() > timeLimit) {
                    solver.cancelSearch();
                }
            });

            solver.PropagationCompleted.subscribe_unhandled([&stopWatch, &solver, &solution, this](auto &&) {
                if ((hardTimeLimit or solution) and stopWatch.elapsed<std::chrono::milliseconds>() > timeLimit) {
                    solver.cancelSearch();
                }
            });
            std::vector<BooleanVar<ST>> variables;
            std::vector<ST> weights;
            variables.reserve(gnnOutput.size());
            weights.reserve(gnnOutput.size());
            for (auto weight: gnnOutput | elements<1>) {
                auto x = solver.newBoolean();
                solver.addToSearch(x);
                variables.emplace_back(x);
                weights.emplace_back(static_cast<ST>(weight * FixPointPrec));
            }

            const auto &cb = proxy.getSolver().clauses;
            if (cb.size() != 0) {
                auto lookup = detail::makeLookup(solver, variables, proxy.getSolver().numLiteral(),
                                                 gnnOutput | elements<0>);
                auto clauses = iota(0ul, cb.size()) | transform([&cb](auto idx) { return cb[idx]; }) |
                               filter([](auto ptr) { return nullptr != ptr; });
                for (const auto c : clauses) {
                    auto clause = *c | transform(lookup) | common;
                    solver.clauses.add(clause.begin(), clause.end());
                }
            }

            solver.post(AtMost(static_cast<ST>(numLiterals), variables));
            auto objective = Sum(variables, weights);
            stopWatch.start();
            solver.maximize(objective);
            if (not solver.boolean.hasSolution()) {
                proxy.fail();
                return 0;
            }

            for (auto [valid, lit]: iterators::zip(
                    variables | transform([&solver](const auto &x) { return solver.boolean.value(x); }),
                    gnnOutput | elements<0>)) {
                if (valid) {
                    cache.emplace_back(lit);
                }
            }

            proxy.makeAssumptions(cache);
            return cache.size();
        }
    };

    /**
     * @brief Variant wrapper for fix policies
     * @tparam Ts
     */
    template<typename ...Ts>
    struct VariantFix : public std::variant<Ts...> {
        DYNAMIC_DISPATCH_VOID(reset, EMPTY)

        DYNAMIC_DISPATCH_FORWARD(select, EMPTY)
    };

    template<concepts::scalar T>
    using FixPolicy = VariantFix<BestN, GreedyFix<T>, OptimalFix<T>>;

    /**
     * @brief GNN based relaxation policy. Uses precedence predictor to fix most probable precedences
     * @tparam T timing type
     * @tparam R resource type
     */
    template<concepts::scalar T, SchedulingResource R>
    class GNNBackbonePredictor {
        GNNPrecedencePredictor<T, R> predictor;
        const Solver<T> &solver;
        mutable tempo::util::Profiler profiler{};
        const PolicyConfig config;
        std::vector<std::pair<Literal<T>, double>> gnnCache;
        FixPolicy<T> fixPolicy{};
        unsigned failCount = 0;
        unsigned solverFailCount = 0;
        double fixRatio;
        std::size_t numFixed = 0;

        auto calcFailRatio(auto fails) const noexcept {
            return static_cast<double>(fails - solverFailCount) / predictor.numLiterals();
        }

        auto decayFactor(auto failRatio) const {
            using enum DecayMode;
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

    public:
        GNNBackbonePredictor(const GNNBackbonePredictor &) = default;
        GNNBackbonePredictor(GNNBackbonePredictor &&) = default;

        /**
         * Dtor. Prints timing information if solver verbosity allows it
         */
        ~GNNBackbonePredictor() {
            if (solver.getOptions().verbosity >= Options::YACKING) {
                profiler.printAll<std::chrono::milliseconds>(std::cout);
            }
        }

        /**
         * Ctor. Performs call to the GNN based on the current solver state
         * @param solver solver in use
         * @param modelLocation path to the model file
         * @param featureExtractorConfigLocation path to the feature extractor config
         * @param problemInstance problem instance
         * @param relaxationRatio percentage of search literals to relay (in [0, 1])
         * @param relaxationDecay decay factor applied to the relaxation ratio on each relaxation fail
         * @param minCertainty minimum GNN prediction certainty
         */
        GNNBackbonePredictor(const Solver<T> &solver, const fs::path &modelLocation,
                             const fs::path &featureExtractorConfigLocation,
                             const SchedulingProblemHelper<T, R> &problemInstance,
                             const PolicyConfig &config) :
                predictor(modelLocation, featureExtractorConfigLocation, problemInstance,
                          problemInstance.getSearchLiterals(solver)),
                          solver(solver), config(config), fixRatio(config.fixRatio) {
            if (config.fixRatio < 0 or config.fixRatio > 1) {
                throw std::runtime_error("invalid fix ratio");
            }

            if (config.decay < 0 or config.decay >= 1) {
                throw std::runtime_error("invalid decay");
            }

            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << config << std::endl;
            }

            using enum AssumptionMode;
            switch(config.assumptionMode) {
                case GreedySkip:
                    fixPolicy.template emplace<GreedyFix<T>>(false);
                    break;
                case GreedyInverse:
                    fixPolicy.template emplace<GreedyFix<T>>(true);
                    break;
                case Optimal:
                    fixPolicy.template emplace<OptimalFix<T>>(100, 50, false);
                    break;
                case AssumptionMode::SingleShot:
                    break;
            }
        }

        /**
         * Repair interface. Fixes a number of literals based on the current fix ratio
         * @tparam AssumptionInterface
         * @param s
         */
        template<heuristics::assumption_interface AssumptionInterface>
        void fix(AssumptionInterface &s) {
            const auto numLits = maxNumLiterals();
            if (numLits == 0) {
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "repair");
            numFixed = fixPolicy.select(s, numLits, failCount, gnnCache);
            if (solver.getOptions().verbosity >= Options::YACKING) {
                if (s.getState() == heuristics::AssumptionState::Fail) {
                    std::cout << "-- failed to fix literals\n";
                } else {
                    std::cout << "-- fixing " << numFixed << " / " << predictor.numLiterals()
                              << " literals" << (failCount > 0 ? " after fail\n" : "\n");
                }
            }
        }

        /**
         * Call this method when the last relaxation was a success. Currently does nothing
         */
        void notifySuccess(unsigned fails) {
            const auto failRatio = calcFailRatio(fails);
            if (failRatio > config.maxFailRatio and config.decreaseOnSuccess) {
                fixRatio *= decayFactor(failRatio);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- decreasing fix ratio to " << fixRatio * 100
                              << "% after too many solver fails" << std::endl;
                }
            } else if (failRatio < config.minFailRatio) {
                fixRatio = std::min(1.0, fixRatio / config.decay);
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- increasing fix ratio to " << fixRatio * 100
                              << "%" << std::endl;
                }
            }

            solverFailCount = fails;
            failCount = 0;
        }

        /**
         * Runs the GNN inference and fills the assumption cache
         */
        void runInference() {
            using namespace std::views;
            if (maxNumLiterals() == 0) {
                gnnCache.clear();
                return;
            }

            tempo::util::ScopeWatch sw(profiler, "gnn lns policy update");
            predictor.updateConfidence(solver);
            auto lits = predictor.getLiterals();
            auto selection = lits | take_while([m = config.minCertainty](auto &tpl) { return std::get<1>(tpl) > m; })
                             | common;
            gnnCache = std::vector(selection.begin(), selection.end());
            if (solver.getOptions().verbosity >= Options::YACKING) {
                std::cout << "-- Updating GNN confidence values" << std::endl;
            }
        }

        /**
         * Maximum number of literals the predictor will try to fix
         * @return upper bound on number of literals to fix
         */
        auto maxNumLiterals() const noexcept {
           return static_cast<std::size_t>(predictor.numLiterals() * fixRatio);
        }

        /**
         * Call this function to reset confidences and caches and rerun inference. This is useful when the gnn will be
         * applied to a new region in the search space
         */
        void notifyNewRegion() {
            failCount = 0;
            predictor.reinitialize(solver);
            fixRatio = config.fixRatio;
            runInference();
            fixPolicy.reset();
        }

        /**
         * Indicates whether a new region should be explored
         * @return true if ratio of literals to fix is lower than exhaustion threshold, false otherwise
         */
        bool exhausted() const noexcept {
            return (std::min(numFixed, maxNumLiterals()) /
                   static_cast<double>(predictor.numLiterals())) < config.exhaustionThreshold;
        }

        /**
         * Call this method when the last relaxation lead to an UNSAT problem. Updates the GNN prediction and decreases
         * the relaxation ratio
         */
        void notifyFailure(unsigned fails) {
            if (++failCount > config.retryLimit) {
                fixRatio *= decayFactor(calcFailRatio(fails));
                if (maxNumLiterals() > 0 and solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << std::setprecision(2) << "-- setting fix ratio = "
                              << fixRatio * 100 << "%" << std::endl;
                }
            } else {
                if (solver.getOptions().verbosity >= Options::YACKING) {
                    std::cout << "-- backbone prediction failed, retrying more carefully\n";
                }
            }

            solverFailCount = fails;
        }
    };
}

#endif //TEMPO_GNNBACKBONEPREDICTOR_HPP
