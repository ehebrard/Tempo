/**
* @author Tim Luchterhand
* @date 05.09.24
* @brief utility program that adds number of tasks and variables to the info file for all problems under a folder
*/

#include <iostream>
#include <filesystem>

#include "data_generation.hpp"
#include "util/scheduling_helpers.hpp"
#include "util/Options.hpp"
#include "helpers/cli.hpp"
#include "helpers/VarImportanceRunner.hpp"


namespace fs = std::filesystem;

int main(int argc, char **argv) {
    using namespace tempo;
    auto opt = cli::parseOptions(argc, argv);
    const fs::path problems(opt.instance_file);
    for (const auto &entry : fs::directory_iterator(problems)) {
        if (not fs::is_directory(entry)) {
            continue;
        }

        std::cout << "processing " << entry.path() << std::endl;
        auto meta = getInfo(entry);
        opt.instance_file = getInstance(entry);
        auto [_, _1, _2, _3, nTasks] = loadSchedulingProblem(opt);
        auto [dp, status] = tempo::loadDataPoint(entry, 0, true);
        if (status != DataPointStatus::Valid) {
            std::cerr << "failed to load problem " << entry.path() << std::endl;
            continue;
        }

        auto numVariables = VarImportanceRunner(dp.problem, opt, dp.solution).getLiterals().size();
        meta["numberOfTasks"] = nTasks;
        meta["numberOfVariables"] = numVariables;
        serialization::serializeToFile(meta, entry.path() / InfoFileName);
    }

}