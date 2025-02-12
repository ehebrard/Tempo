/**
* @author Tim Luchterhand
* @date 24.09.24
* @brief
*/

#include <filesystem>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <Iterators.hpp>

#include "data_generation.hpp"
#include "nn/torch_types.hpp"
#include "nn/tensor_utils.hpp"
#include "nn/GraphBuilder.hpp"
#include "util/Matrix.hpp"
#include "util/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"
#include "../helpers/git_sha.hpp"

auto createGraphBuilder(const std::string &featConfig, const tempo::Options &options) {
    auto [_, problem, _1, _2, _3, _4] = loadSchedulingProblem(options);
    return tempo::nn::GraphBuilder(featConfig, std::move(problem));
}


int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;

    std::string saveTo;
    std::string featConfig;
    auto options = cli::parseOptions(argc, argv,
                                     cli::ArgSpec("save-to", "Where to save the features", true, saveTo),
                                     cli::ArgSpec("feat-config", "Path to the feature extractor config", true,
                                                  featConfig));
    if (not fs::is_directory(options.instance_file)) {
        std::cerr << "specify the path to the directory with the data points" << std::endl;
        std::exit(1);
    }

    const auto mainDir = fs::path(options.instance_file);
    options.instance_file = getInstance(options.instance_file);
    const auto inputsDir = fs::path(saveTo) / InputsName;
    const auto labelDir = fs::path(saveTo) / OutputsName;
    fs::create_directories(inputsDir);
    fs::create_directory(labelDir);
    auto solutions = getSolutions(mainDir);
    if (solutions.empty()) {
        std::cerr << "no solutions found" << std::endl;
        std::exit(1);
    }

    const auto optimalSol = solutions.rbegin()->second;
    nlohmann::json optSolPayload;
    {
        auto [s, p, _, _1, _2, _3] = loadSchedulingProblem(options);
        try {
            s->set(leq(p.schedule().duration.id(), optimalSol.objective));
            loadBranch(*s, optimalSol.decisions);
        } catch (const Failure<int> &) {
            std::cerr << "solution " << optimalSol.id << " of problem '" << mainDir.filename() << "' is inconsistent"
                      << std::endl;
            std::exit(1);
        }
        auto taskDistances = p.getTaskDistances(*s);
        Matrix<int> network(p.tasks().size(), p.tasks().size(), std::move(taskDistances));
        optSolPayload["taskNetwork"] = std::move(network);
        optSolPayload["objective"] = optimalSol.objective;
        optSolPayload["id"] = optimalSol.id;
    }

    auto builder = createGraphBuilder(featConfig, options);
    fs::copy(featConfig, saveTo, fs::copy_options::overwrite_existing);
    std::vector<nn::InputGraph> graphs;
    graphs.reserve(solutions.size());
    for (const auto &[id, sol] : solutions) {
        auto [solver, problem, _, _1, _2, _3] = loadSchedulingProblem(options);
        try {
            solver->set(leq(problem.schedule().duration.id(), sol.objective));
            loadBranch(*solver, sol.decisions);
        } catch (const Failure<int> &) {
            std::cerr << "solution " << id << " of problem '" << mainDir.filename() << "' is inconsistent" << std::endl;
            continue;
        }

        graphs.emplace_back(builder.getGraph(nn::makeSolverState(problem.getTaskDistances(*solver), *solver)));
    }

    for (auto [id, graph] : iterators::const_enumerate(graphs)) {
        const auto destination = inputsDir / (GraphName + std::to_string(id));
        fs::create_directories(destination);
        nn::util::saveGraph(graph, destination);
        const auto lDir = labelDir / (LabelName + std::to_string(id));
        fs::create_directories(lDir);
        serialization::serializeToFile(optSolPayload, lDir / LabelFileName);
    }

    nlohmann::json meta;
    meta["bestSolutionKnown"] = optimalSol.objective;
    meta["numSubProblems"] = graphs.size();
    meta["numRootProblems"] = 0;
    meta["date"] = shell::getTimeStamp();
    meta["commit"] = GitSha;

    serialization::serializeToFile(meta, fs::path(saveTo) / InfoFileName);
    return 0;
}