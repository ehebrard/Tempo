/**
* @author Tim Luchterhand
* @date 08.08.24
* @brief
*/

#include <ranges>
#include <Iterators.hpp>

#include "VarImportanceRunner.hpp"
#include "scheduling_helpers.hpp"

VarImportanceRunner::VarImportanceRunner(ser::PartialProblem partialProblem, tempo::Options options,
                                         const ser::Solution<Time> &optSol) : problem(std::move(partialProblem)),
                                                                              options(std::move(options)),
                                                                              optimum(optSol.objective) {
    auto problemInstance = loadSchedulingProblem(this->options);
    auto &s = *problemInstance.solver;
    const auto &p = problemInstance.instance;
    try {
        loadBranch(s, problem.decisions);
    } catch (const tempo::Failure<Time> &) {
        inconsistent = true;
        return;
    }

    tempo::var_t maxVar = 0;
    for (auto var : s.getBranch()) {
        maxVar = std::max(maxVar, var);
        if (s.boolean.hasSemantic(var)) {
            auto edge = s.boolean.getEdge(true, var);
            if (p.hasVariable(edge.from) and p.hasVariable(edge.to)) {
                searchLiterals.emplace_back(s.boolean.getLiteral(true, var));
            }
        }
    }

    literalCache.resize(maxVar * 2 + 2, false);
    for (const auto &[var, val] : optSol.decisions) {
        if (val and var <= maxVar) {
            auto lit = s.boolean.getLiteral(val, var);
            literalCache.at(lit) = true;
        }
    }
}

auto VarImportanceRunner::getLiterals() const noexcept -> const std::vector<tempo::Literal<Time>> & {
    return searchLiterals;
}

double VarImportanceRunner::averageSearchEffort() const noexcept {
    if (searchLiterals.empty()) {
        return 0;
    }

    return averageNumberOfDecisions() / static_cast<double>(searchLiterals.size());
}

double VarImportanceRunner::averageNumberOfDecisions() const noexcept {
    if (totalNumDecisions == 0 and numSearches == 0) {
        return 0;
    }

    return static_cast<double>(totalNumDecisions) / numSearches;
}

auto VarImportanceRunner::run(tempo::Literal<Time> lit) -> Result {
    using enum SchedulerState;
    if (isInconsistent()) {
        throw std::runtime_error("inconsistent sub problem");
    }

    if (literalCache.at(lit)) {
        return {Valid, optimum};
    }

    auto problemInstance = loadSchedulingProblem(options);
    auto &s = *problemInstance.solver;
    const auto &p = problemInstance.instance;
    loadBranch(s, problem.decisions);
    if (s.boolean.satisfied(lit) or s.boolean.falsified(lit)) {
        return {AlreadyDecided, {}};
    }

    try {
        s.set(lit);
    } catch(const tempo::Failure<Time> &) {
        return {Unsat, {}};
    }

    s.minimize(p.schedule().duration);
    totalNumDecisions += s.num_choicepoints;
    ++numSearches;
    if (s.boolean.hasSolution() and s.numeric.hasSolution()) {
        const auto makeSpan = s.numeric.lower(p.schedule().duration);
        if (makeSpan == optimum) {
            for (auto [litId, cacheVal] : iterators::enumerate(literalCache)) {
                cacheVal = cacheVal or s.boolean.bestSolution().at(litId);
            }
        }

        return {Valid, makeSpan};
    }

    return {Unsat, {}};
}

bool VarImportanceRunner::isInconsistent() const noexcept {
    return inconsistent;
}
