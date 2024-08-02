/**
* @author Tim Luchterhand
* @date 30.07.24
* @brief Program for verification of data points
*/

#include <filesystem>
#include <iostream>
#include <optional>

#include "util/serialization.hpp"
#include "data_generation.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "helpers/cli.hpp"

auto setDifference(const tempo::serialization::Branch &a, const tempo::serialization::Branch &b) {
    tempo::serialization::Branch out;
    std::ranges::set_difference(a, b, std::back_inserter(out));
    return out;
}

auto findMatchingSolution(const tempo::serialization::PartialProblem &p,
                          const std::filesystem::path &mainDir)
                          -> std::pair<std::optional<tempo::serialization::Solution<int>>,
                          tempo::serialization::Branch> {
    using namespace tempo;
    const auto solutions = getSolutions(mainDir);
    auto decisions = p.decisions;
    std::ranges::sort(decisions);
    auto diff = setDifference(decisions, solutions.at(p.associatedSolution).decisions);
    if (diff.empty()) {
        return {solutions.at(p.associatedSolution), std::move(diff)};
    }

    for (const auto &sol : solutions) {
        if (setDifference(decisions, sol.decisions).empty()) {
            return {sol, std::move(diff)};
        }
    }

    return {{}, std::move(diff)};
}


int main(int argc, char **argv) {
    using namespace tempo;
    namespace ser = tempo::serialization;
    namespace fs = std::filesystem;
    unsigned dpId;

    auto options = cli::parseOptions(argc, argv, cli::ArgSpec("dp-id", "id of the data point", true, dpId));
    const auto mainDir = fs::path(options.instance_file);
    std::optional<ser::PartialProblem> partialProblem;
    for (const auto &file : fs::directory_iterator(mainDir / Serializer<>::SubProblemDir)) {
        if (file.path().filename() ==
            Serializer<>::SubProblemBaseName + std::to_string(dpId) + Serializer<>::FileExtension) {
            partialProblem = ser::deserializeFromFile<ser::PartialProblem>(file);
            break;
        }
    }

    if (not partialProblem.has_value()) {
        std::cerr << "sub problem with id " << dpId << " not found" << std::endl;
        std::exit(1);
    }

    ser::Solution<int> solution;
    try {
        solution = getSolutions(mainDir).at(partialProblem->associatedSolution);
    } catch (const std::out_of_range&) {
        std::cerr << "associated solution " << partialProblem->associatedSolution << " not found" << std::endl;
        std::exit(1);
    }

    auto [matching, diff] = findMatchingSolution(*partialProblem, mainDir);
    if (not diff.empty()) {
        std::cout << "indicated solution is no superset of partial problem. Conflicts: ";
        for (auto [var, val]: diff) {
            std::cout << "(" << var << ", " << val << "); ";
        }

        std::cout << std::endl;
    }

    if (matching.has_value()) {
        std::cout << "found solution that is super set of constraints in partial problem. ID " << matching->id
                  << " vs indicated " << partialProblem->associatedSolution << std::endl;
    } else {
        std::cout << "none of the solutions is a super set of the constraints in the partial problems" << std::endl;
    }

    options.instance_file = getInstance(options.instance_file);
    auto [solver, problem, _] = loadSchedulingProblem(options);
    try {
        loadBranch(*solver, partialProblem->decisions);
    } catch(const Failure<int> &f) {
        std::cout << "inconsistent snapshot: " << f.what() << std::endl;
        std::exit(2);
    }

    solver->minimize(problem.schedule().duration);
    if (not solver->numeric.hasSolution()) {
        std::cout << "no solution for data point" << std::endl;
        std::exit(3);
    }
    const int makeSpan = solver->numeric.lower(problem.schedule().duration);
    if (makeSpan != solution.objective) {
        std::cout << "locally optimal solution differs from expected solutions: makespan " << makeSpan << " vs "
                  << solution.objective << std::endl;
        std::exit(4);
    } else {
        std::cout << "found solution with expected makespan" << std::endl;
    }

    return 0;
}