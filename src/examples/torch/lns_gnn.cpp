/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief
*/


#include <iostream>
#include <string>

#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"
#include "../helpers/git_sha.hpp"
#include "nn/GNNRelaxationPolicy.hpp"
#include "util/Profiler.hpp"

int main(int argc, char **argv) {
    using namespace tempo;
    std::string gnnLocation;
    std::string featureExtractorConf;
    auto confidenceThresh = 0.0;
    auto relaxationRatio = 0.0;
    auto relaxationDecay = 0.0;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", false, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", false,
                                              featureExtractorConf),
                                 cli::ArgSpec("confidence", "minimum confidence of GNN", false, confidenceThresh),
                                 cli::ArgSpec("ratio", "percentage of literals to relax", false, relaxationRatio),
                                 cli::ArgSpec("decay", "relaxation ratio decay on failure", false, relaxationDecay));
    auto [solver, problem, optSol, _] = loadSchedulingProblem(opt);
    nn::GNNRelaxationPolicy policy(*solver, gnnLocation, featureExtractorConf, problem, relaxationRatio,
                                   relaxationDecay, confidenceThresh);

    MinimizationObjective objective(problem.schedule().duration);
    util::StopWatch sw;
    solver->largeNeighborhoodSearch(objective, policy);
    auto [start, end] = sw.getTiming();
    if (solver->numeric.hasSolution()) {
        auto makespan = solver->numeric.lower(problem.schedule().duration);
        std::cout << "-- makespan " << makespan << std::endl;
        if (optSol.has_value() and makespan > *optSol) {
            std::cout << "-- suboptimal solution! Optimum is " << *optSol << std::endl;
        }
    }

    std::cout << "-- total duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}