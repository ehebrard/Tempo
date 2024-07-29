/**
* @author Tim Luchterhand
* @date 24.07.24
* @brief program for generating datapoints for nn training
*/

#include <string>

#include "helpers/cli.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "data_generation.hpp"


int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    namespace fs = std::filesystem;

    std::string saveTo;
    auto options = cli::parseOptions(argc, argv,
                                     cli::ArgSpec("save-to", "Where to save the data points", true, saveTo));
    auto [solver, problem] = loadSchedulingProblem(options);
    auto schedule = problem.schedule();
    auto heuristic = make_compound_heuristic(make_variable_heuristic(*solver), TightestSolutionGuided(0, 0));
    solver->setBranchingHeuristic(std::move(heuristic));
    const auto problemName = fs::path(options.instance_file).filename();
    const auto destinationFolder = fs::path(saveTo) / problemName;
    Tracer tracer(*solver, schedule, destinationFolder);
    solver->minimize(schedule.end);
    return 0;
}


