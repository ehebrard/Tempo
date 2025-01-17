/**
* @author Tim Luchterhand
* @date 17.01.25
* @file calc_backbone.cpp
* @brief Calculates the backbone of a set of solutions (decisions that never change). More precisely, calculates the
* variable bias. A bias of -1 or 1 indicates that the variable is part of the backbone
*/

#include <algorithm>
#include <ranges>
#include <vector>
#include <string>

#include "util/Lookup.hpp"
#include "util/serialization.hpp"
#include "data_generation.hpp"
#include "helpers/cli.hpp"
#include "helpers/git_sha.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "helpers/shell.hpp"


#define JSONIFY(JSON, DATA) JSON[#DATA] = DATA;

int main(int argc, char **argv) {
    using namespace tempo;
    namespace v = std::views;
    double qualityThreshold = 0;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("quality", "Quality threshold [0, ∞)", true, qualityThreshold));
    const auto mainDir = std::move(opt.instance_file);
    opt.instance_file = getInstance(mainDir);
    const auto solutions = getSolutions(mainDir);
    if (solutions.empty()) {
        std::cerr << "No solutions found" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const auto optimum = solutions.rbegin()->second.objective;
    auto selectionIt = std::ranges::find_if(solutions, [optimum, qualityThreshold](const auto &pair) {
        auto quality = static_cast<double>(pair.second.objective) / optimum - 1;
        return quality <= qualityThreshold;
    });
    if (selectionIt == solutions.end()) {
        std::cerr << std::setprecision(2) << "No solutions found that are at most " << qualityThreshold * 100 <<
                "% worse as the optimum" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const auto selection = std::ranges::subrange(selectionIt, solutions.end()) | v::elements<1>;
    const auto selectionSize = std::ranges::distance(selection);
    Lookup variableBias(solutions.rbegin()->second.decisions | v::elements<0>, 0.0);
    for (const auto &sol : selection) {
        for (auto [var, val] : sol.decisions) {
            variableBias[var] += 2 * val - 1;
        }
    }

    unsigned backBoneSize = 0;
    unsigned unimportantSize = 0;
    const auto solutionSize = solutions.rbegin()->second.decisions.size();
    std::vector<std::string> taskEdgeBias; std::vector<std::pair<var_t, double>> varBias;
    taskEdgeBias.reserve(solutionSize); varBias.reserve(solutionSize);
    const EdgeMapper mapper(opt);
    for (auto var : solutions.rbegin()->second.decisions | v::elements<0>) {
        auto bias = variableBias[var] / static_cast<double>(selectionSize);
        backBoneSize += std::abs(bias) == 1;
        unimportantSize += bias == 0;
        varBias.emplace_back(var, bias);
        auto [f, t] = mapper.getTaskEdge(mapper.getSolver().boolean.getLiteral(true, var));
        taskEdgeBias.emplace_back(std::to_string(f) + " <-> " + std::to_string(t) + " = " + std::to_string(bias));
    }

    nlohmann::json json;
    JSONIFY(json, taskEdgeBias)
    JSONIFY(json, varBias)
    JSONIFY(json, unimportantSize)
    JSONIFY(json, backBoneSize)
    JSONIFY(json, qualityThreshold)
    JSONIFY(json, selectionSize)
    json["date"] = shell::getTimeStamp();
    json["commit"] = GitSha;
    std::cout << json.dump(__JSON_INDENT__) << std::endl;
    return 0;
}
