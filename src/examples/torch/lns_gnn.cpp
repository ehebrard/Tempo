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
#include "heuristics/DRPolicy.hpp"
#include "heuristics/relaxation_policy_factories.hpp"

namespace h = tempo::heuristics;

using RP = h::RelaxationPolicy<Time, ResourceConstraint>;


int main(int argc, char **argv) {
    using namespace tempo;
    std::string gnnLocation;
    std::string featureExtractorConf;
    nn::PolicyConfig config;
    h::RelaxationPolicyParams destroyParameters{.relaxRatio = 0.9, .ratioDecay = 1, .numScheduleSlices = 4};
    h::RelaxPolicy destroyType;
    double sporadicIncrement = 0.001;
    double exhaustionProbability = 0.1;
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
                                 cli::ArgSpec("exhaustion-threshold",
                                              "fix rate lower threshold at which a new region should be chosen",
                                              false, config.exhaustionThreshold),
                                 cli::ArgSpec("exhaustion-prob",
                                              "probability of choosing a new region even when GNN is not exhausted",
                                              false, exhaustionProbability),
                                 cli::ArgSpec("destroy-ratio", "percentage of literals to relax", false,
                                              destroyParameters.relaxRatio),
                                 cli::ArgSpec("destroy-decay", "decay applied to the fix ratio (inverse destroy ratio)", false,
                                              destroyParameters.ratioDecay),
                                 cli::ArgSpec("num-slices", "number of schedule slices", false,
                                              destroyParameters.numScheduleSlices),
                                 cli::ArgSpec("sporadic-increment", "probability increment on fail for root search",
                                              false, sporadicIncrement),
                                 cli::ArgSpec("fix-ratio", "percentage of literals to relax", false,
                                              config.fixRatio),
                                 cli::ArgSpec("decay", "relaxation ratio reactivity on failure", false,
                                              config.decay),
                                 cli::ArgSpec("assumption-mode", "how to make assumptions", false,
                                              config.assumptionMode),
                                 cli::SwitchSpec("decrease-on-success", "whether to decrease fix rate even on success",
                                                 config.decreaseOnSuccess, false),
                                 cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio",
                                              false, config.retryLimit),
                                 cli::ArgSpec("decay-mode", "relaxation ratio decay mode on failure", false,
                                              config.decayMode),
                                 cli::ArgSpec("destroy-mode", "destroy policy type", true, destroyType),
                                 cli::ArgSpec("threads", "GNN inference threads", false,
                                              numThreads));
    auto problemInfo = loadSchedulingProblem(opt);
    torch::set_num_threads(numThreads);
    nn::GNNBackbonePredictor gnnRepair(*problemInfo.solver, gnnLocation, featureExtractorConf, problemInfo.instance,
                                       config);

    std::cout << "-- root search probability increment " << sporadicIncrement << std::endl;
    std::cout << "-- exhaustion probability " << exhaustionProbability << std::endl;
    h::GenericDestroyPolicy<Time, RP> destroy(
            h::make_relaxation_policy(destroyType, problemInfo.instance.tasks(), problemInfo.constraints,
                                      destroyParameters));
    std::cout << "-- using destroy policy " << destroyType << std::endl;
    auto policy = heuristics::make_sporadic_root_search(sporadicIncrement,
                                                        heuristics::make_RD_policy(destroy, gnnRepair,
                                                                                   exhaustionProbability));
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
