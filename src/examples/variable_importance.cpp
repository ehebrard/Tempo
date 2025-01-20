/**
* @author Tim Luchterhand
* @date 06.08.24
* @brief Calculates the importance of binary decisions for all search variables in a given (sub) problem
*/

#include <string>
#include <limits>
#include <unordered_map>
#include <nlohmann/json.hpp>

#include "helpers/shell.hpp"
#include "helpers/cli.hpp"
#include "helpers/git_sha.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "data_generation.hpp"
#include "util/KillHandler.hpp"
#include "helpers/VarImportanceRunner.hpp"

#define JSONIFY(JSON, DATA) JSON[#DATA] = DATA;

namespace fs = std::filesystem;


int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    int subId = -1;
    bool root;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("sub-number", "id of the sub problem", false, subId),
                                 cli::SwitchSpec("root", "calculate edge importance for a root instance instead", root,
                                                 false));
    const auto mainDir = opt.instance_file;
    opt.instance_file = getInstance(mainDir);
    DataPoint dataPoint;
    if (root) {
        auto solutions = getSolutions(mainDir);
        if (solutions.empty()) {
            std::cerr << "No solutions found" << std::endl;
            std::exit(1);
        }

        dataPoint.solution = std::move(solutions.rbegin()->second);
    } else {
        DataPointStatus status;
        std::tie(dataPoint, status) = loadDataPoint(mainDir, subId, false);
        if (status == DataPointStatus::ProblemNotFound) {
            std::cerr << "sub problem with id " << subId << " could not be found" << std::endl;
            std::exit(1);
        } else if (status == DataPointStatus::SolutionNotFound) {
            std::cerr << "associated solution" << subId << " could not be found" << std::endl;
            std::exit(1);
        }
    }

    VarImportanceRunner runner(dataPoint.problem, opt, dataPoint.solution);
    if (runner.isInconsistent()) {
        std::cout << "snapshot is inconsistent" << std::endl;
        std::exit(2);
    }

    unsigned numIndifferent = 0;
    unsigned numDecisions = 0;
    unsigned numDecided = 0;
    unsigned numEdges = 0;
    const auto optimum = dataPoint.solution.objective;
    std::unordered_map<var_t, double> varImportance;
    const EdgeMapper mapper(opt);
    std::vector<std::string> taskEdgeImportance;
    for (auto lit : runner.getLiterals()) {
        if (varImportance.contains(lit.variable())) {
            //This should never happen
            std::cout << "duplicate literal in search variables" << std::endl;
            std::exit(2);
        }

        using enum SchedulerState;
        auto posResult = runner.run(lit);
        if (KillHandler::instance().signalReceived()) {
            break;
        }

        ++numEdges;
        if (posResult.state == AlreadyDecided) {
            ++numDecided;
            continue;
        }

        auto negResult = runner.run(~lit);
        if (negResult.state == AlreadyDecided) {
            ++numDecided;
            continue;
        }

        ++numDecisions;
        if (posResult.result.value_or(0) != optimum and negResult.result.value_or(0) != optimum) {
            std::cout << "ERROR: " << lit << " (variable " << lit.variable() << ")"
                      << " is a choice point but the optimal solution cannot be found from this state" << std::endl;
            std::exit(2);
        }

        if (posResult.state == Unsat and negResult.state == Unsat) {
            std::cout << "ERROR: " << lit << " (variable " << lit.variable() << ")"
                      << "cannot be set! Inconsistent datapoint" << std::endl;
            std::exit(2);
        }

        constexpr auto Inf = std::numeric_limits<double>::infinity();
        const auto resPos = static_cast<std::optional<double>>(posResult.result).value_or(Inf);
        const auto resNeg = static_cast<std::optional<double>>(negResult.result).value_or(Inf);
        const auto larger = std::max(resPos, resNeg);
        const auto smaller = std::min(resPos, resNeg);
        const auto quality = 1 - static_cast<double>(smaller) / larger;
        numIndifferent += quality == 0;
        varImportance[lit.variable()] = quality;
        auto [f, t] = mapper.getTaskEdge(lit);
        taskEdgeImportance.emplace_back(
                std::to_string(f) + " <-> " + std::to_string(t) + "=" + std::to_string(quality));
    }

    auto unimportant = static_cast<double>(numIndifferent) / numDecisions;
    auto decided = static_cast<double>(numDecided) / numEdges;
    if (std::isnan(unimportant)) { unimportant = 0; }
    if (std::isnan(decided)) { decided = 1; }
    nlohmann::json json;
    JSONIFY(json, varImportance);
    JSONIFY(json, taskEdgeImportance);
    JSONIFY(json, unimportant);
    JSONIFY(json, decided);
    json["problem"] = fs::path(mainDir).filename();
    json["searchEffort"] = runner.averageSearchEffort();
    json["avgDecisions"] = runner.averageNumberOfDecisions();
    json["date"] = shell::getTimeStamp();
    json["commit"] = GitSha;
    std::cout << json.dump(__JSON_INDENT__) << std::endl;
    return 0;
}
