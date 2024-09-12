/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#include <iostream>
#include <string>
#include <vector>

#include "nn/GNNPrecedencePredictor.hpp"
#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"


void setLiterals(tempo::Solver<int> &solver, const std::vector<tempo::Literal<int>> &literals) {
    if (solver.getOptions().verbosity >= tempo::Options::NORMAL) {
        std::cout << "-- setting " << literals.size() << " literals\n";
    }

    for (auto lit : literals) {
        auto branch = solver.getBranch();
        if (not branch.safe_has(lit.variable())) {
            continue;
        }

        if (solver.getOptions().verbosity >= tempo::Options::NORMAL) {
            std::cout << lit << " ";
        }
        solver.set(lit);
    }

    std::cout << std::endl;
}


int main(int argc, char **argv) {
    using namespace tempo;
    std::string gnnLocation;
    std::string featureExtractorConf;
    unsigned numberOfIterations = 0;
    double confidenceThresh = 1;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", true, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", true,
                                              featureExtractorConf),
                                 cli::ArgSpec("confidence", "minimum confidence of GNN", true, confidenceThresh),
                                 cli::ArgSpec("iterations", "number of precedence predictor runs", false,
                                              numberOfIterations));

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

    nn::GNNPrecedencePredictor precedencePredictor(gnnLocation, featureExtractorConf, problem, std::move(literals));
    if (numberOfIterations > 0) {
        precedencePredictor.updateConfidence(*solver);
    }

    setLiterals(*solver, precedencePredictor.getLiterals(confidenceThresh));
    solver->minimize(schedule.duration);
    if (solver->numeric.hasSolution()) {
        auto makespan = solver->numeric.lower(schedule.duration);
        std::cout << "-- makespan " << makespan << std::endl;
        if (optSol.has_value() and makespan > *optSol) {
            std::cout << "-- suboptimal solution! Optimum is " << *optSol << std::endl;
        }
    }

    return 0;
}