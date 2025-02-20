/**
* @author Tim Luchterhand
* @date 20.02.25
* @file lns_region_to_taillard.cpp
* @brief Converts LNS regions to taillard format:
* #jobs #machines #precedences ub
* <job configs>
* ...
* <newline>
* machine-t1 t1 machine-t2 t2
* ...
*/

#include <filesystem>
#include <string>
#include <ranges>
#include <Iterators.hpp>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "util/SchedulingProblemHelper.hpp"
#include "util/serialization.hpp"
#include "heuristics/LNS/relaxation_evaluators.hpp"
#include "helpers/cli.hpp"

namespace fs = std::filesystem;
namespace ser = tempo::serialization;
using TaskPrecedence = std::pair<unsigned, unsigned>;

auto getTaskResourceMap(const tempo::SchedulingProblemHelper<tempo::Time, tempo::Resource> &instance) -> std::vector<unsigned> {
    std::vector<unsigned> taskToResource(instance.tasks().size());
    for (unsigned task = 0; task < instance.tasks().size(); ++task) {
        for (auto [idx, r] : iterators::enumerate(instance.resources())) {
            if (std::ranges::find(r.tasks(), task) != r.tasks().end()) {
                taskToResource[task] = idx;
                break;
            }
        }
    }

    return taskToResource;
}

auto splitInstance(const fs::path &instanceFile) -> std::pair<std::string, std::string> {
    std::ifstream ifs(instanceFile);
    if (not ifs.is_open()) {
        std::cerr << "Could not open instance file " << instanceFile << " for reading." << std::endl;
        std::exit(1);
    }

    std::string header;
    if (not std::getline(ifs, header)) {
        std::cerr << "Could not read header from instance file " << instanceFile << std::endl;
        std::exit(1);
    }

    std::ostringstream oss;
    oss << ifs.rdbuf();
    return {std::move(header), oss.str()};
}

int main(int argc, char **argv) {
    using namespace tempo;
    using namespace std::views;
    bool onlySat;
    std::string regionFiles;
    std::string destination;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("regions", "path to folder file with regions", true, regionFiles),
                                 cli::ArgSpec("out", "path to folder where sub problems are stored", true, destination),
                                 cli::SwitchSpec("sat", "whether to only include SAT instances", onlySat, false));
    if (not fs::is_directory(regionFiles)) {
        std::cerr << "The region directory " << destination << " does not exist." << std::endl;
        std::exit(1);
    }
    if (not fs::is_directory(destination)) {
        std::cerr << "The destination directory " << destination << " does not exist." << std::endl;
        std::exit(1);
    }

    const auto instanceName = fs::path(opt.instance_file).filename().string();
    const auto regionFileName = fs::path(regionFiles) / instanceName;
    if (not fs::exists(regionFileName)) {
        std::cerr << "The region file " << regionFileName << " could not be found in " << regionFiles << std::endl;
        std::exit(1);
    }

    const auto regions = ser::deserializeFromFile<std::vector<lns::Region<Time>>>(regionFileName);
    const auto problem = loadSchedulingProblem(opt);
    const auto &instance = problem.instance;

    std::vector<std::vector<TaskPrecedence>> taskPrecedences;
    auto filteredRegions = regions | filter([onlySat](const auto &r) { return not onlySat or r.SAT; });
    for (const auto &region: filteredRegions) {
        std::vector<TaskPrecedence> tp;
        const auto &v2t = instance.getMapping();
        for (const auto &dc: region.assumptions | filter([&instance](const auto &dc) {
            return instance.hasVariable(dc.from) and instance.hasVariable(dc.to);
        })) {
            tp.emplace_back(v2t(dc.to), v2t(dc.from));
        }

        taskPrecedences.emplace_back(std::move(tp));
    }

    const auto [header, rest] = splitInstance(opt.instance_file);
    const auto t2m = getTaskResourceMap(instance);
    for (auto [idx, taskPrecs, region] : iterators::zip_enumerate(taskPrecedences, filteredRegions)) {
        const auto name = instanceName + "_" + std::to_string(idx) + "_" + (region.SAT ? "sat" : "unsat");
        const auto file = fs::path(destination) / name;
        std::ofstream ofs(file);
        if (not ofs.is_open()) {
            std::cerr << "Could not open file " << file << " for writing." << std::endl;
            std::exit(1);
        }

        ofs << header << " " << taskPrecs.size() << " " << region.ub << "\n" << rest << "\n";
        for (auto [f, t]: taskPrecs) {
            ofs << t2m.at(f) << " " << f << " " << t2m.at(t) << " " << t << "\n";
        }
    }

    return 0;
}
