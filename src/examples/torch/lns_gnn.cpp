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
                                 cli::ArgSpec("update-threshold", "failure rate threshold at which to update the GNN",
                                              false, config.confidenceUpdateFailRatio),
                                 cli::ArgSpec("ratio", "percentage of literals to relax", false,
                                              config.fixRatio),
                                 cli::ArgSpec("decay", "relaxation ratio reactivity on failure", false,
                                              config.decay),
                                 cli::SwitchSpec("careful", "whether to make careful assumptions after failure",
                                                  config.carefulAssumptions, false),
                                 cli::SwitchSpec("decrease-on-success", "whether to decrease fix rate even on success",
                                                 config.decreaseOnSuccess, false),
                                 cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio",
                                              false, config.retryLimit),
                                 cli::ArgSpec("decay-mode", "relaxation ratio decay mode on failure", false,
                                              config.decayMode),
                                 cli::ArgSpec("threads", "GNN inference threads", false,
                                              numThreads));
    auto problemInfo = loadSchedulingProblem(opt);
    torch::set_num_threads(numThreads);
    nn::GNNBackbonePredictor policy(*problemInfo.solver, gnnLocation, featureExtractorConf, problemInfo.instance,
                                    config);
    MinimizationObjective objective(problemInfo.instance.schedule().duration);
    util::StopWatch sw;
    problemInfo.solver->largeNeighborhoodSearch(objective, policy);
    auto [start, end] = sw.getTiming();
    if (problemInfo.solver->numeric.hasSolution()) {
        auto makespan = problemInfo.solver->numeric.lower(problemInfo.instance.schedule().duration);
        std::cout << "-- makespan " << makespan << std::endl;
        if (problemInfo.optimalSolution.has_value() and makespan > *problemInfo.optimalSolution) {
            std::cout << "-- suboptimal solution! Optimum is " << *problemInfo.optimalSolution << std::endl;
        }
    }

    std::cout << "-- total duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}
