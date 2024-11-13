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
    lns::RelaxationPolicyParams policyParams;
    lns::RelaxPolicy policyType;
    auto options = cli::parseOptions(argc, argv,
                                     cli::ArgSpec("save-to", "Where to save the data points", true, saveTo),
                                     cli::SwitchSpec("lns", "Whether to use LNS", useLNS, false),
                                     cli::ArgSpec("relax-decay", "relaxation ratio decay",
                                                  false, policyParams.ratioDecay, 0.5),
                                     cli::ArgSpec("relax-ratio", "initial relaxation ratio",
                                                  false, policyParams.relaxRatio, 0.5),
                                     cli::ArgSpec("relax-slices", "number of schedule slices",
                                                  false, policyParams.numScheduleSlices, 4),
                                     cli::ArgSpec("lns-policy", "lns relaxation policy", false, policyType,
                                                  lns::RelaxPolicy::RandomTasks));
    auto [solver, problem, constraints, opt, nTasks] = loadSchedulingProblem(options);
    auto schedule = problem.schedule();
    auto heuristic = make_compound_heuristic(make_variable_heuristic(*solver), TightestSolutionGuided(0, 0));
    solver->setBranchingHeuristic(std::move(heuristic));
    const auto problemName = fs::path(options.instance_file).filename();
    const auto destinationFolder = fs::path(saveTo) / problemName;
    DataGenerator dataGenerator(*solver, schedule, destinationFolder);
    if (useLNS) {
        MinimizationObjective<int> objective(schedule.duration);
        auto policy = lns::make_relaxation_policy(policyType, problem.tasks(), constraints, policyParams);
        std::cout << "-- using relaxation policy " << policyType << std::endl;
        solver->largeNeighborhoodSearch(objective, policy);
    } else {
        solver->minimize(schedule.duration);
    }

    std::optional<int> makespan;
    if (solver->numeric.hasSolution()) {
        makespan = solver->numeric.lower(schedule.duration);
    }

    fs::copy(options.instance_file, destinationFolder / ProblemFileName, fs::copy_options::overwrite_existing);
    nlohmann::json meta;
    meta["date"] = shell::getTimeStamp();
    meta["commit"] = GitSha;
    meta["numSubProblems"] = dataGenerator.problemCount();
    meta["numSolutions"] = dataGenerator.solutionCount();
    meta["bestSolutionKnown"] = optMin(opt, makespan);
    meta["numberOfTasks"] = nTasks;
    serialization::serializeToFile(meta, destinationFolder / InfoFileName);
    return 0;
}


