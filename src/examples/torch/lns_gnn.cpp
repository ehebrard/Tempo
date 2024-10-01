/**
* @author Tim Luchterhand
* @date 18.09.24
* @brief
*/


#include <iostream>
#include <string>
#include <thread>

#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"
#include "../helpers/git_sha.hpp"
#include "nn/GNNBackbonePredictor.hpp"
#include "util/Profiler.hpp"

int main(int argc, char **argv) {
    using namespace tempo;
    std::string gnnLocation;
    std::string featureExtractorConf;
    nn::PolicyConfig config;
    unsigned numThreads = std::max(1u, std::thread::hardware_concurrency() / 2);
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", false, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", false,
                                              featureExtractorConf),
                                 cli::ArgSpec("confidence", "minimum confidence of GNN", false,
                                              config.minCertainty),
                                 cli::ArgSpec("min-fail", "lower bound solver failure rate", false,
                                              config.minFailRatio),
                                 cli::ArgSpec("max-fail", "upper bound solver failure rate", false,
                                              config.maxFailRatio),
                                 cli::ArgSpec("ratio", "percentage of literals to relax", false,
                                              config.relaxationRatio),
                                 cli::ArgSpec("decay", "relaxation ratio decay on failure", false,
                                              config.relaxationDecay),
                                 cli::SwitchSpec("careful", "whether to make careful assumptions after failure",
                                                  config.carefulAssumptions, false),
                                 cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio", false,
                                              config.retryLimit),
                                 cli::ArgSpec("threads", "GNN inference threads", false,
                                              numThreads));
    auto [solver, problem, optSol, _] = loadSchedulingProblem(opt);
    torch::set_num_threads(numThreads);
    nn::GNNBackbonePredictor policy(*solver, gnnLocation, featureExtractorConf, problem, config);
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
