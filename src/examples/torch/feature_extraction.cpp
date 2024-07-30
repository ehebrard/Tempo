/**
* @author Tim Luchterhand
* @date 29.07.24
* @brief program used for extracting torch features for GNN training
*/

#include <filesystem>
#include <vector>
#include <string>
#include <ranges>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <Iterators.hpp>

#include "data_generation.hpp"
#include "nn/torch_types.hpp"
#include "nn/tensor_utils.hpp"
#include "nn/GraphBuilder.hpp"
#include "util/KillHandler.hpp"
#include "util/Matrix.hpp"
#include "../helpers/scheduling_helpers.hpp"
#include "../helpers/cli.hpp"

namespace fs = std::filesystem;
constexpr auto RootName = "root";
constexpr auto GraphName = "graph";
constexpr auto InputsName = "inputs";
constexpr auto OutputsName = "outputs";
constexpr auto LabelFileName = "task_network";
constexpr auto LabelName = "label";
constexpr auto InfoFileName = "info.json";


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
                                     cli::SwitchSpec("filter-duplicates", "Whether to filter duplicate graphs",
                                                     filterDuplicates, true));
    if (not fs::is_directory(options.instance_file)) {
        std::cerr << "specify the path to the directory with the data points" << std::endl;
        std::exit(1);
    }


    const auto mainDir = fs::path(options.instance_file);
    const auto problemsDir = mainDir / Serializer<>::SubProblemDir;
    options.instance_file = getInstance(options.instance_file);
    const auto inputsDir = fs::path(saveTo) / InputsName;
    const auto labelDir = fs::path(saveTo) / OutputsName;
    const auto rootInputDir = inputsDir / RootName;
    fs::create_directories(rootInputDir);
    fs::create_directory(labelDir);

    auto [solver, problem, opt] = loadSchedulingProblem(options);
    nn::GraphBuilder builder(featConfig, problem);
    fs::copy(featConfig, saveTo, fs::copy_options::overwrite_existing);

    // create root instance features
    solver->propagate();
    nn::util::saveGraph(builder.getGraph(nn::makeSolverState(problem.getTaskDistances(*solver), *solver)),
                        rootInputDir);

    // Serialize solutions
    auto solutions = getSolutions(mainDir);
    std::vector<nlohmann::json> solutionsPayloads;
    for (const auto &sol : solutions) {
        auto [s, p, _] = loadSchedulingProblem(options);
        loadBranch(*s, sol.decisions);
        auto taskDistances = p.getTaskDistances(*s);
        Matrix<int> network(p.tasks().size(), p.tasks().size(), taskDistances);
        nlohmann::json j;
        j["taskNetwork"] = network;
        j["objective"] = sol.objective;
        j["id"] = sol.id;
        solutionsPayloads.emplace_back(std::move(j));
    }

    // create input features
    std::vector<std::tuple<nn::InputGraph, unsigned>> problemGraphs;
    unsigned numDuplicates = 0;
    for (const auto &file : fs::directory_iterator(problemsDir)) {
        using std::views::elements;
        if (not file.is_regular_file()) {
            continue;
        }

        auto partial = serialization::deserializeFromFile<serialization::PartialProblem>(file);
        auto [s, p, _] = loadSchedulingProblem(options);
        loadBranch(*s, partial.decisions);
        auto newGraph = builder.getGraph(nn::makeSolverState(p.getTaskDistances(*s), *s));
        bool equal = false;
        if (filterDuplicates) {
            for (const auto &g : problemGraphs | elements<0>) {
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
            problemGraphs.emplace_back(std::move(newGraph), partial.associatedSolution);
        }
    }

    std::cout << "filtered out " << numDuplicates << " duplicates" << std::endl;
    for (auto [id, g] : iterators::const_enumerate(problemGraphs)) {
        const auto &[graph, solId] = g;
        const auto destination = inputsDir / (GraphName + std::to_string(id));
        fs::create_directories(destination);
        nn::util::saveGraph(graph, destination);
        const auto &solPayload = solutionsPayloads.at(solId);
        if (solPayload.at("id").get<unsigned>() != solId) {
            std::cerr << "this should not happen. Maybe there are solution files missing?" << std::endl;
            std::exit(-1);
        }

        const auto lDir = labelDir / (LabelName + std::to_string(id));
        fs::create_directories(lDir);
        serialization::serializeToFile(solPayload, lDir / LabelFileName);
    }

    for (const auto &solPayload : solutionsPayloads) {
        const auto destination = labelDir / (RootName + std::to_string(solPayload.at("id").get<unsigned>()));
        fs::create_directories(destination);
        serialization::serializeToFile(solPayload, destination / LabelFileName);
    }

    auto bestSolution = opt.value_or(solutions.back().objective);
    nlohmann::json meta;
    meta["bestSolutionKnown"] = bestSolution;
    meta["numSubProblems"] = problemGraphs.size();
    meta["numRootProblems"] = solutions.size();
    serialization::serializeToFile(meta, fs::path(saveTo) / InfoFileName);
    return 0;
}