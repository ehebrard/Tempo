/**
* @author Tim Luchterhand
* @date 10.07.24
* @brief
*/

#include "heuristics/GNNValueHeuristics.hpp"
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
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("gnn-loc", "Location of the GNN model", true, gnnLocation),
                                 cli::ArgSpec("feat-config", "Location of the feature extractor config", true,
                                              featureExtractorConf));

    auto [solver, problem, _] = loadSchedulingProblem(opt);
    auto schedule = problem.schedule();
    using PType = decltype(problem);
    using GNNH = util::ProfiledHeuristic<GNNFullGuidance<PType::Time, PType::Resource>>;
    util::Profiler profiler;
    GNNH valBranching(profiler, opt.polarity_epsilon, gnnLocation, featureExtractorConf, std::move(problem));
    auto varBranching = make_variable_heuristic(*solver);
    solver->setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(valBranching)));
    solver->minimize(schedule.duration);
    if (solver->numeric.hasSolution()) {
        std::cout << "-- makespan " << solver->numeric.lower(schedule.duration) << std::endl;
    }

    profiler.printAll<std::chrono::milliseconds>(std::cout);
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}