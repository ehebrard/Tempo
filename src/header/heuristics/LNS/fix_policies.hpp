/**
* @author Tim Luchterhand
* @date 07.11.24
* @brief Contains literal selection policies for relax / repair LNS-policies
*/

#ifndef TEMPO_FIX_POLICIES_HPP
#define TEMPO_FIX_POLICIES_HPP

#include <vector>
#include <ranges>

#include "util/traits.hpp"
#include "Literal.hpp"
#include "util/Options.hpp"
#include "Model.hpp"
#include "util/Profiler.hpp"
#include "relaxation_interface.hpp"
#include "relaxation_policies.hpp"
#include "Solver.hpp"

namespace tempo::lns {

    PENUM(AssumptionMode, BestN, GreedySkip, GreedyInverse, Optimal)

    /**
     * @brief Literal ordering type. Literals are ordered by their weight either in ascending or descending order or
     * they are unordered
     */
    enum class OrderType {
        Unordered,
        Ascending,
        Descending
    };

    /**
     * @brief Fix policy that simply takes the best N edges predicted by the GNN
     * @tparam OT Assumes that literals are ordered by their weights according to the given ordering type
     */
    template<OrderType OT>
    struct BestN {

        constexpr void reset() noexcept {};

        template<concepts::scalar T, assumption_interface AI, concepts::scalar C>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned,
                           const std::vector<std::pair<Literal<T>, C>> & weightedLiterals) {
            using namespace std::views;
            switch (OT) {
                case OrderType::Unordered: {
                    auto copy = weightedLiterals;
                    std::ranges::sort(copy, {}, [](const auto &elem) { return -std::get<1>(elem); });
                    proxy.makeAssumptions(copy | elements<0> | take(numLiterals));
                    break;
                }
                case OrderType::Ascending:
                    proxy.makeAssumptions(weightedLiterals | reverse | take(numLiterals) | elements<0>);
                    break;
                case OrderType::Descending:
                    proxy.makeAssumptions(weightedLiterals | take(numLiterals) | elements<0>);
                    break;
            }
            return std::min(numLiterals, weightedLiterals.size());
        }
    };

    /**
     * @brief Fix policy that tries to fix N edges in a greedy fashion
     * @tparam T timing type
     * @tparam OT Assumes that literals are ordered by their weights according to the given ordering type
     */
    template<concepts::scalar T, OrderType OT>
    class GreedyFix {
        std::vector<Literal<T>> cache;
        bool inverse;

        auto getLiterals(const auto &weightedLits) const {
            using namespace std::views;
            using enum OrderType;
            if constexpr (OT == Unordered) {
                auto copy = weightedLits;
                std::ranges::sort(copy, {}, [](const auto &elem) { return -std::get<1>(elem); });
                return std::ranges::owning_view(std::move(copy)) | elements<0>;
            } else if constexpr (OT == Ascending) {
                return weightedLits | reverse | elements<0>;
            } else {
                return weightedLits | elements<0>;
            }
        }
    public:
        /**
         * Ctor
         * @param inverse whether to set the inverse literal on fail
         */
        explicit constexpr GreedyFix(bool inverse) noexcept: inverse(inverse) {}

        void reset() noexcept {
            cache.clear();
        }

        template<assumption_interface AI, concepts::scalar C>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned numFails,
                           const std::vector<std::pair<Literal<T>, C>> & weightedLiterals) {
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
            for (auto lit: getLiterals(weightedLiterals)) {
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

    template<assumption_interface AI, concepts::scalar C>
    std::size_t select(AI &proxy, std::size_t numLiterals, unsigned numFails,
                       const std::vector<std::pair<Literal<T>, C>> & weightedLiterals) {
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
        variables.reserve(weightedLiterals.size());
        weights.reserve(weightedLiterals.size());
        for (auto weight: weightedLiterals | elements<1>) {
            auto x = solver.newBoolean();
            solver.addToSearch(x);
            variables.emplace_back(x);
            weights.emplace_back(static_cast<ST>(weight * FixPointPrec));
        }

        const auto &cb = proxy.getSolver().clauses;
        if (cb.size() != 0) {
            auto lookup = detail::makeLookup(solver, variables, proxy.getSolver().numLiteral(),
                                             weightedLiterals | elements<0>);
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
                weightedLiterals | elements<0>)) {
            if (valid) {
                cache.emplace_back(lit);
            }
        }

        proxy.makeAssumptions(cache);
        return cache.size();
    }
};

    template<concepts::scalar T>
    class TaskFix {
        static_assert(traits::always_false_v<T>, "Implementation incomplete, do not use!");
        std::vector<std::pair<Interval<T>, double>> tasks;
        lns::detail::TaskVarMap<T> map;
        std::vector<Literal<T>> cache{};
        std::size_t varsPerTask = 0;
        std::size_t lastNumTasks = 0;

        void calcWeights(const std::vector<std::pair<Literal<T>, double>> & weightedLiterals) {
            for (auto &[t, confidence]: tasks) {
                confidence = 0;
                for (const auto &var : map(t)) {
                    auto res = std::ranges::find_if(weightedLiterals, var, [&var](const auto &pair) {
                        return pair.first.variable() == var.id();
                    });

                    if (res != weightedLiterals.end()) {
                        confidence += res->second;
                    }
                }
            }

            std::ranges::sort(tasks, {}, [](const auto &pair) { return -pair.second; });
        }
    public:

        template<concepts::typed_range<Interval<T>> Tasks,  resource_range RR>
        TaskFix(const Tasks &tasks, const RR &resources): map(tasks, resources) {
            this->tasks.reserve(std::ranges::size(tasks));
            for (const auto &t : tasks) {
                varsPerTask += map(t).size();
                this->tasks.emplace_back(t, 0.0);
            }
        }

        void reset() noexcept {
            cache.clear();
        }

        template<assumption_interface AI>
        std::size_t select(AI &proxy, std::size_t numLiterals, unsigned numFails,
                           const std::vector<std::pair<Literal<T>, double>> & weightedLiterals) {
            const auto numTasks = numLiterals / varsPerTask;
            if (numTasks == 0) {
                return 0;
            }

            if (lastNumTasks == numTasks and not cache.empty()) {
                proxy.makeAssumptions(cache);
                return cache.size();
            }

            if (cache.empty()) {
                calcWeights(weightedLiterals);
            }
        }
    };


    /**
     * @brief Variant wrapper for fix policies
     * @tparam Ts
     */
    template<typename ...Ts>
    struct VariantFix : public std::variant<Ts...> {
        using std::variant<Ts...>::variant;

        DYNAMIC_DISPATCH_VOID(reset, EMPTY)

        DYNAMIC_DISPATCH_FORWARD(select, EMPTY)
    };
}

#endif //TEMPO_FIX_POLICIES_HPP
