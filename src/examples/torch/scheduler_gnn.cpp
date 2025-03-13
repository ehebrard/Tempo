/**
* @author Tim Luchterhand
* @date 10.07.24
* @brief
*/

#include "nn/GNNValueHeuristics.hpp"
#include "heuristics/heuristic_factories.hpp"
#include "util/parsing/jsp.hpp"
#include "util/Profiler.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/shell.hpp"
#include "../helpers/git_sha.hpp"


int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    std::string gnnLocation;
    std::string featureExtractorConf;
    nn::DispatcherConfig dispatcherConfig {
        .failIncrement = 0.0001, .successDecrement = 0.0005, .restartIncrement = 0.05, .solutionDecrement = 0.05,
        .maxFillRate = 0.2, .heatIncrement = 0.5, .heatDecay = 0.95, .heatLowerThreshold = 0.1
    };
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", true, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", true,
                                              featureExtractorConf),
                                 cli::ArgSpec("gnn-fail-increment", "dispatcher fail increment in %", false,
                                              dispatcherConfig.failIncrement),
                                 cli::ArgSpec("gnn-success-decrement", "dispatcher success decrement in %", false,
                                              dispatcherConfig.successDecrement),
                                 cli::ArgSpec("gnn-restart-increment", "dispatcher restart increment in %", false,
                                              dispatcherConfig.restartIncrement),
                                 cli::ArgSpec("gnn-solution-decrement", "dispatcher solution decrement in %", false,
                                              dispatcherConfig.solutionDecrement),
                                 cli::ArgSpec("gnn-max-fill-rate", "dispatcher max fill rate in %", false,
                                              dispatcherConfig.maxFillRate),
                                 cli::ArgSpec("gnn-heat-increment", "dispatcher heat increment on inference in %",
                                              false, dispatcherConfig.heatIncrement),
                                 cli::ArgSpec("gnn-heat-decay", "dispatcher heat decay", false,
                                              dispatcherConfig.heatDecay));

    auto [solver, problem, _, _1, _2, _3] = loadSchedulingProblem(opt);
    auto schedule = problem.schedule();
    using PType = decltype(problem);
    using GNNH = util::ProfiledHeuristic<nn::GNNFullGuidance<PType::Time, PType::Resource>>;
    util::Profiler profiler;
    GNNH valBranching(profiler, opt.polarity_epsilon, gnnLocation, featureExtractorConf, *solver, dispatcherConfig,
                      std::move(problem));
    auto varBranching = make_variable_heuristic(*solver);
    solver->setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(valBranching)));
    solver->minimize(schedule.duration);
    if (solver->numeric.hasSolution()) {
        std::cout << "-- makespan " << solver->numeric.solutionLower(schedule.duration) << std::endl;
    }

    profiler.printAll<std::chrono::milliseconds>(std::cout);
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}