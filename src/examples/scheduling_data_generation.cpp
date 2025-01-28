/**
* @author Tim Luchterhand
* @date 24.07.24
* @brief program for generating datapoints for nn training
*/

#include <string>
#include <nlohmann/json.hpp>

#include "helpers/cli.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "helpers/shell.hpp"
#include "data_generation.hpp"
#include "helpers/git_sha.hpp"
#include "heuristics/LNS/relaxation_policy_factories.hpp"

template<typename T>
auto optMin(const std::optional<T> &lhs, const std::optional<T> &rhs) -> std::optional<T> {
    if (not lhs.has_value() and not rhs.has_value()) { return lhs; }
    if (lhs.has_value() and rhs.has_value()) { return std::min(*lhs, *rhs); }
    return lhs.has_value() ? lhs : rhs;
}

int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    namespace fs = std::filesystem;

    std::string saveTo;
    bool useLNS;
    lns::RelaxationPolicyParams policyParams{
        .decayConfig = lns::PolicyDecayConfig(), .numScheduleSlices = 4, .allTaskEdges = false
    };
    lns::RelaxPolicy policyType;
    double sporadicIncrement = 0.001;
    auto options = cli::parseOptions(argc, argv,
                                     cli::ArgSpec("save-to", "Where to save the data points", true, saveTo),
                                     cli::SwitchSpec("lns", "Whether to use LNS", useLNS, false),
                                     cli::ArgSpec("fix-decay", "relaxation ratio decay",
                                                  false, policyParams.decayConfig.decay, 0.5),
                                     cli::ArgSpec("fix-ratio", "initial relaxation ratio",
                                                  false, policyParams.decayConfig.fixRatio, 0.5),
                                     cli::ArgSpec("decay-min-fail",
                                                  "lower bound solver failure rate for ratio decay config",
                                                  false, policyParams.decayConfig.minFailRatio),
                                     cli::ArgSpec("decay-max-fail",
                                                  "upper bound solver failure rate for ratio decay config",
                                                  false, policyParams.decayConfig.maxFailRatio),
                                     cli::SwitchSpec("decay-on-success", "whether to decrease fix rate even on success",
                                                     policyParams.decayConfig.decreaseOnSuccess, false),
                                     cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio",
                                                  false, policyParams.decayConfig.retryLimit),
                                     cli::ArgSpec("decay-mode", "relaxation ratio decay mode on failure", false,
                                                  policyParams.decayConfig.decayMode),
                                     cli::ArgSpec("relax-slices", "number of schedule slices",
                                                  false, policyParams.numScheduleSlices, 4),
                                     cli::ArgSpec("lns-policy", "lns relaxation policy", false, policyType,
                                                  lns::RelaxPolicy::RandomTasks),
                                     cli::ArgSpec("sporadic-increment", "probability increment on fail for root search",
                                                  false, sporadicIncrement),
                                     cli::SwitchSpec("fix-all-task-edges",
                                                     "whether to fix all task edges or only those between fixed tasks",
                                                     policyParams.allTaskEdges, false)
    );
    auto problemInfo = loadSchedulingProblem(options);
    const auto schedule = problemInfo.instance.schedule();
    auto heuristic = make_compound_heuristic(make_variable_heuristic(*problemInfo.solver), TightestSolutionGuided(0, 0));
    problemInfo.solver->setBranchingHeuristic(std::move(heuristic));
    const auto problemName = fs::path(options.instance_file).filename();
    const auto destinationFolder = fs::path(saveTo) / problemName;
    DataGenerator dataGenerator(*problemInfo.solver, schedule, destinationFolder);
    if (useLNS) {
        MinimizationObjective<int> objective(schedule.duration);
        auto policy = lns::make_sporadic_root_search(sporadicIncrement, lns::make_relaxation_policy(
                                                         policyType, problemInfo.instance.tasks(),
                                                         problemInfo.constraints, policyParams));
        std::cout << "-- using relaxation policy " << policyType << std::endl;
        problemInfo.solver->largeNeighborhoodSearch(objective, policy);
    } else {
        problemInfo.solver->minimize(schedule.duration);
    }

    std::optional<int> makespan;
    if (problemInfo.solver->numeric.hasSolution()) {
        makespan = problemInfo.solver->numeric.lower(schedule.duration);
    }

    fs::copy(options.instance_file, destinationFolder / ProblemFileName, fs::copy_options::overwrite_existing);
    nlohmann::json meta;
    meta["date"] = shell::getTimeStamp();
    meta["commit"] = GitSha;
    meta["numSubProblems"] = dataGenerator.problemCount();
    meta["numSolutions"] = dataGenerator.solutionCount();
    meta["bestSolutionKnown"] = optMin(problemInfo.optimalSolution, makespan);
    meta["numberOfTasks"] = problemInfo.numTasks;
    serialization::serializeToFile(meta, destinationFolder / InfoFileName);
    return 0;
}


