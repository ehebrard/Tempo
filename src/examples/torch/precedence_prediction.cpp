/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#include <iostream>
#include <string>
#include <variant>
#include <vector>
#include <optional>
#include <ranges>

#include "nn/GNNPrecedencePredictor.hpp"
#include "heuristics/TightestPrecedencePredictor.hpp"
#include "util/Profiler.hpp"
#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"


template<tempo::concepts::ttyped_range<tempo::Literal> L>
void setLiterals(tempo::Solver<int> &solver, L &&literals) {
    std::size_t numLits = 0;
    for (auto lit : std::forward<L>(literals)) {
        auto branch = solver.getBranch();
        if (not branch.safe_has(lit.variable())) {
            continue;
        }

        if (solver.getOptions().verbosity >= tempo::Options::NORMAL) {
            std::cout << lit << " ";
        }
        solver.set(lit);
        ++numLits;
    }

    std::cout << std::endl;
    if (solver.getOptions().verbosity >= tempo::Options::NORMAL) {
        std::cout << "-- setting " << numLits << " literals" << std::endl;
    }

}

enum class PredictorType {
    GNN, Tightest
};

using GNNPredictor = tempo::nn::GNNPrecedencePredictor<int, DisjunctiveResource<int>>;
using TightestPredictor = tempo::heuristics::TightestPrecedencePredictor<int>;

using Predictor = std::variant<TightestPredictor, GNNPredictor>;

int main(int argc, char **argv) {
    using namespace tempo;
    std::string gnnLocation;
    std::string featureExtractorConf;
    int pType = 0;
    unsigned numberOfIterations = 0;
    auto confidenceThresh = 0.0;
    unsigned numLiterals = 0;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", false, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", false,
                                              featureExtractorConf),
                                 cli::ArgSpec("confidence", "minimum confidence of GNN", false, confidenceThresh),
                                 cli::ArgSpec("num-lits", "number of literals to set", false, numLiterals),
                                 cli::ArgSpec("iterations", "number of precedence predictor runs", false,
                                              numberOfIterations),
                                 cli::ArgSpec("predictor-type", "predictor type to use", false, pType));
    const PredictorType predictorType{pType};
    auto [solver, problem, optSol, _] = loadSchedulingProblem(opt);
    auto schedule = problem.schedule();
    std::vector<Literal<int>> literals;
    for (auto var : solver->getBranch()) {
        if (solver->boolean.hasSemantic(var)) {
            auto edge = solver->boolean.getEdge(true, var);
            if (problem.hasVariable(edge.from) and problem.hasVariable(edge.to)) {
                literals.emplace_back(solver->boolean.getLiteral(true, var));
            }
        }
    }

    std::optional<Predictor> predictor;
    if (predictorType == PredictorType::GNN) {
        predictor.emplace(std::in_place_type<GNNPredictor>, gnnLocation, featureExtractorConf, problem, std::move(literals));
    } else if (predictorType == PredictorType::Tightest) {
        predictor.emplace(std::in_place_type<TightestPredictor>, std::move(literals));
    } else {
        std::cerr << "unknown predictor type" << std::endl;
        std::exit(1);
    }

    util::StopWatch sw;
    if (numberOfIterations > 0 and numLiterals > 0) {
        auto o = opt;
        o.verbosity = Options::SILENT;
        auto [s, p, _1, _2] = loadSchedulingProblem(o);
        s->setBranchingHeuristic(heuristics::make_compound_heuristic(heuristics::RandomVariableSelection{},
                                                                     heuristics::RandomBinaryValue{}));
        s->PropagationCompleted.subscribe_unhandled(
                [i = 0u, numberOfIterations, &predictor](auto &state) mutable {
            if (++i > numberOfIterations) {
                KillHandler::instance().kill();
            } else {
                std::visit([&state](auto &pred) {pred.updateConfidence(state);}, *predictor);
            }
        });

        s->minimize(p.schedule().duration);
        KillHandler::instance().reset();

        using namespace std::views;
        auto allLits = std::visit([](const auto &pred) { return pred.getLiterals(); }, *predictor);
        setLiterals(*solver,
                    allLits |
                    filter([confidenceThresh](const auto &tpl) { return std::get<1>(tpl) > confidenceThresh; }) |
                    take(numLiterals) | elements<0>);
    }

    solver->minimize(schedule.duration);
    auto [start, end] = sw.getTiming();
    if (solver->numeric.hasSolution()) {
        auto makespan = solver->numeric.lower(schedule.duration);
        std::cout << "-- makespan " << makespan << std::endl;
        if (optSol.has_value() and makespan > *optSol) {
            std::cout << "-- suboptimal solution! Optimum is " << *optSol << std::endl;
        }
    }

    std::cout << "total duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms";
    return 0;
}