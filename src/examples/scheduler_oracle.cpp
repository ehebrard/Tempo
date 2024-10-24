/**
* @author Tim Luchterhand
* @date 13.08.24
* @brief runs the scheduler using the perfect oracle value heuristic
*/

#include "helpers/scheduling_helpers.hpp"
#include "heuristics/PerfectValueOracle.hpp"
#include "helpers/cli.hpp"
#include "helpers/shell.hpp"
#include "helpers/git_sha.hpp"
#include "data_generation.hpp"
#include "util/Profiler.hpp"


int main(int argc, char **argv) {
    using namespace tempo;
    namespace ser = tempo::serialization;
    namespace fs = std::filesystem;
    using namespace heuristics;

    unsigned solutionId;
    auto options = cli::parseOptions(argc, argv, cli::ArgSpec("solution-id", "id of the solution", true, solutionId));
    const auto mainDir = fs::path(options.instance_file);
    options.instance_file = getInstance(mainDir);
    const auto [dataPoint, status] = loadDataPoint(mainDir, solutionId, true);
    const auto &[_, solution] = dataPoint;
    if (status == DataPointStatus::SolutionNotFound) {
        std::cerr << "associated solution " << solutionId << " not found" << std::endl;
        std::exit(1);
    }

    auto [solver, problem, _1, _2, _3] = loadSchedulingProblem(options);
    using Oracle = util::ProfiledHeuristic<PerfectValueHeuristic>;
    util::Profiler profiler;
    Oracle valueOracle(profiler, options.polarity_epsilon, solution);
    solver->setBranchingHeuristic(make_compound_heuristic(make_variable_heuristic(*solver),
                                                          std::move(valueOracle)));
    solver->minimize(problem.schedule().duration);
    if (solver->numeric.hasSolution()) {
        std::cout << "-- makespan " << solver->numeric.lower(problem.schedule().duration) << std::endl;
    }

    profiler.printAll<std::chrono::microseconds>(std::cout);
    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}