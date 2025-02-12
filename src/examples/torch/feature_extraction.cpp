/**
* @author Tim Luchterhand
* @date 29.07.24
* @brief program used for extracting torch features for GNN training
*/

#include <filesystem>
#include <vector>
#include <string>
#include <ranges>
#include <nlohmann/json.hpp>
#include <Iterators.hpp>

#include "data_generation.hpp"
#include "nn/torch_types.hpp"
#include "nn/tensor_utils.hpp"
#include "nn/GraphBuilder.hpp"
#include "util/KillHandler.hpp"
#include "util/Matrix.hpp"
#include "util/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"
#include "../helpers/shell.hpp"
#include "../helpers/git_sha.hpp"

namespace fs = std::filesystem;

template<tempo::concepts::scalar T>
bool distanceSanityCheck(const tempo::Matrix<T> &distances) {
    for (std::size_t i = 0; i < distances.numRows(); ++i) {
        for (std::size_t j = 0; j < distances.numColumns(); ++j) {
            if (i == j) {
                continue;
            }

            auto dist = distances(i, j) + distances(j, i);
            if (dist < 0) {
                return false;
            }
        }
    }

    return true;
}

int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;

    std::string saveTo;
    std::string featConfig;
    bool filterDuplicates;
    auto options = cli::parseOptions(argc, argv,
                                     cli::ArgSpec("save-to", "Where to save the features", true, saveTo),
                                     cli::ArgSpec("feat-config", "Path to the feature extractor config", true,
                                                  featConfig),
                                     cli::SwitchSpec("no-filter-duplicates", "disable duplicate graph filtering",
                                                     filterDuplicates, true));
    if (not fs::is_directory(options.instance_file)) {
        std::cerr << "specify the path to the directory with the data points" << std::endl;
        std::exit(1);
    }


    const auto mainDir = fs::path(options.instance_file);
    options.instance_file = getInstance(options.instance_file);
    const auto inputsDir = fs::path(saveTo) / InputsName;
    const auto labelDir = fs::path(saveTo) / OutputsName;
    const auto rootInputDir = inputsDir / RootName;
    fs::create_directories(rootInputDir);
    fs::create_directory(labelDir);

    Problem problemInfo = loadSchedulingProblem(options);
    nn::GraphBuilder builder(featConfig, problemInfo.instance);
    fs::copy(featConfig, saveTo, fs::copy_options::overwrite_existing);

    // create root instance features
    problemInfo.solver->propagate();
    nn::util::saveGraph(builder.getGraph(
                                nn::makeSolverState(problemInfo.instance.getTaskDistances(*problemInfo.solver), *problemInfo.solver)),
                        rootInputDir);

    // Serialize solutions
    auto solutions = getSolutions(mainDir);
    std::vector<nlohmann::json> solutionsPayloads;
    for (const auto &[id, sol]: solutions) {
        auto [s, p, _, _1, _2, _3] = loadSchedulingProblem(options);
        s->set(leq(p.schedule().duration.id(), sol.objective));
        loadBranch(*s, sol.decisions);
        auto taskDistances = p.getTaskDistances(*s);
        Matrix<int> network(p.tasks().size(), p.tasks().size(), taskDistances);
        if (not distanceSanityCheck(network)) {
            std::cerr << "solution " << id << " of problem '" << mainDir << "' has negative cycle" << std::endl;
            std::exit(-1);
        }

        nlohmann::json j;
        j["taskNetwork"] = network;
        j["objective"] = sol.objective;
        j["id"] = sol.id;
        solutionsPayloads.emplace_back(std::move(j));
    }

    // create input features
    std::vector<std::tuple<nn::InputGraph, unsigned, unsigned>> problemGraphs;
    unsigned numDuplicates = 0;
    for (const auto &[subId, subProblem] : getProblems(mainDir)) {
        using std::views::elements;
        auto [s, p, _, _1, _2, _3] = loadSchedulingProblem(options);
        const auto &sol = solutions.at(subProblem.associatedSolution);
        s->set(leq(p.schedule().duration.id(), sol.objective));
        loadBranch(*s, subProblem.decisions);
        auto newGraph = builder.getGraph(nn::makeSolverState(p.getTaskDistances(*s), *s));
        bool equal = false;
        if (filterDuplicates) {
            for (const auto &g: problemGraphs | elements<0>) {
                if (KillHandler::instance().signalReceived() or nn::util::graphsEquivalent(g, newGraph)) {
                    equal = true;
                    break;
                }
            }
        }

        if (KillHandler::instance().signalReceived()) {
            break;
        }

        numDuplicates += equal;
        if (not equal) {
            problemGraphs.emplace_back(std::move(newGraph), subProblem.associatedSolution, subId);
        }
    }

    std::cout << "filtered out " << numDuplicates << " duplicates" << std::endl;
    for (auto [id, g]: iterators::const_enumerate(problemGraphs)) {
        const auto &[graph, solId, subId] = g;
        const auto destination = inputsDir / (GraphName + std::to_string(id));
        fs::create_directories(destination);
        nlohmann::json ref;
        ref["subId"] = subId;
        ref["solutionId"] = solId;
        nn::util::saveGraph(graph, destination);
        serialization::serializeToFile(ref, destination / GraphReferenceFile);
        const auto &solPayload = solutionsPayloads.at(solId);
        if (solPayload.at("id").get<unsigned>() != solId) {
            std::cerr << "this should not happen. Maybe there are solution files missing?" << std::endl;
            std::exit(-1);
        }

        const auto lDir = labelDir / (LabelName + std::to_string(id));
        fs::create_directories(lDir);
        serialization::serializeToFile(solPayload, lDir / LabelFileName);
    }

    for (const auto &solPayload: solutionsPayloads) {
        const auto destination = labelDir / (RootName + std::to_string(solPayload.at("id").get<unsigned>()));
        fs::create_directories(destination);
        serialization::serializeToFile(solPayload, destination / LabelFileName);
    }

    auto bestSolution = problemInfo.optimalSolution.value_or(solutions.rbegin()->second.objective);
    nlohmann::json meta;
    meta["bestSolutionKnown"] = bestSolution;
    meta["numSubProblems"] = problemGraphs.size();
    meta["numRootProblems"] = solutions.size();
    meta["date"] = shell::getTimeStamp();
    meta["commit"] = GitSha;

    serialization::serializeToFile(meta, fs::path(saveTo) / InfoFileName);
    return 0;
}