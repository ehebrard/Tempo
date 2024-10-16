/**
* @author Tim Luchterhand
* @date 07.10.24
* @brief Use this example to replay a previously recorded LNS policy
*/


#include <iostream>
#include <string>
#include <thread>

#include "helpers/scheduling_helpers.hpp"
#include "helpers/cli.hpp"
#include "helpers/shell.hpp"
#include "helpers/git_sha.hpp"
#include "heuristics/RelaxationPolicy.hpp"

int main(int argc, char **argv) {
    using namespace tempo;
    std::string policyTrace;
    auto opt = cli::parseOptions(argc, argv, cli::ArgSpec("trace", "Location of the policy trace", false, policyTrace));
    auto [solver, problem, _, optSol, _1] = loadSchedulingProblem(opt);
    heuristics::PolicyReplay<int> policy(policyTrace);
    MinimizationObjective objective(problem.schedule().duration);
    solver->largeNeighborhoodSearch(objective, policy);
    if (solver->numeric.hasSolution()) {
        auto makespan = solver->numeric.lower(problem.schedule().duration);
        std::cout << "-- makespan " << makespan << std::endl;
        if (optSol.has_value() and makespan > *optSol) {
            std::cout << "-- suboptimal solution! Optimum is " << *optSol << std::endl;
        }
    }

    std::cout << "-- date: " << shell::getTimeStamp() << std::endl;
    std::cout << "-- commit: " << GitSha << std::endl;
    return 0;
}
