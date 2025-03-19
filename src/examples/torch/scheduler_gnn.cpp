/**
* @author Tim Luchterhand
* @date 10.07.24
* @brief
*/

#include "../helpers/branch_logging.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/git_sha.hpp"
#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/shell.hpp"
#include "heuristics/SingleDecentValueHeuristic.hpp"
#include "heuristics/SolutionGuided.hpp"
#include "heuristics/heuristic_factories.hpp"
#include "nn/GNNDispatcher.hpp"
#include "nn/GNNValueHeuristics.hpp"
#include "util/Profiler.hpp"

int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    std::string gnnLocation;
    std::string featureExtractorConf;
    bool useSolutionGuided = false;
    nn::Dispatch dispatchType = nn::Dispatch::SingleShot;
    std::string recordFile;
    auto opt = cli::parseOptions(
        argc, argv,
        cli::ArgSpec("gnn-loc", "Location of the GNN model", true, gnnLocation),
        cli::ArgSpec("feat-config", "Location of the feature extractor config",
                     true, featureExtractorConf),
        cli::SwitchSpec("solution-guided",
                        "Whether to use solution guided search",
                        useSolutionGuided, false),
        cli::ArgSpec("dispatcher", "dispatcher type", false, dispatchType),
        cli::ArgSpec("record", "record file", false, recordFile, ""));

    if (not recordFile.empty()) {
      opt.restart_policy = "no";
      opt.primal_boost = false;
    }
    std::cout << opt << std::endl;
    auto [solver, problem, _, _1, _2, _3] = loadSchedulingProblem(opt);
    auto schedule = problem.schedule();
    util::Profiler profiler;
    auto varBranching = make_variable_heuristic(*solver);
    auto gnn = make_profiled_heuristic(profiler, nn::make_gnn_value_heuristic(
                                           opt.polarity_epsilon, gnnLocation, featureExtractorConf,
                                           nn::make_dispatcher(dispatchType, *solver), std::move(problem)));
    (*gnn)->runInference(*solver);
    auto valBranching = make_single_decent_heuristic(TightestValue(*solver), std::move(gnn));
    if (useSolutionGuided) {
        std::cout << "-- using solution guided search with GNN" << std::endl;
        auto sg = make_solution_guided_heuristic(*solver, std::move(valBranching));
        solver->setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(sg)));
    } else {
        solver->setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(valBranching)));
    }
    if (not recordFile.empty()) {
      solve(*solver, schedule, recordFile);
    } else {
      solver->minimize(schedule.duration);
    }
    if (solver->numeric.hasSolution()) {
        std::cout << "-- makespan " << solver->numeric.solutionLower(schedule.duration) << std::endl;
    }

    profiler.printAll<std::chrono::milliseconds>(std::cout);
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}