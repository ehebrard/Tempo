/**
* @author Tim Luchterhand
* @date 30.07.24
* @brief
*/

#include "data_generation.hpp"

namespace tempo {

    auto getInstance(const fs::path &problemDir) -> fs::path {
        for (const auto &file : fs::directory_iterator(problemDir)) {
            if (file.is_regular_file() and file.path().extension() == ".txt") {
                return file.path();
            }
        }

        throw std::runtime_error("no problem definition file found in '" + problemDir.string() + "'");
    }

    auto getProblems(const fs::path &problemDir) -> std::vector<serialization::PartialProblem> {
        using namespace tempo::serialization;
        std::vector<PartialProblem> ret;
        const auto dir = problemDir / Serializer<>::SubProblemDir;
        if (not fs::is_directory(dir)) {
            throw std::runtime_error("no problems directory under " + problemDir.string());
        }

        for (const auto &file : fs::directory_iterator(problemDir / Serializer<>::SubProblemDir)) {
            if (file.is_regular_file() and
                file.path().filename().string().starts_with(Serializer<>::SubProblemBaseName)) {
                ret.emplace_back(deserializeFromFile<PartialProblem>(file));
            }
        }

        return ret;
    }
}