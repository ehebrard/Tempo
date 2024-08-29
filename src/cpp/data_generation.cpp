/**
* @author Tim Luchterhand
* @date 30.07.24
* @brief
*/

#include <regex>

#include "data_generation.hpp"

namespace tempo {

    auto getInstance(const fs::path &problemDir) -> fs::path {
        for (const auto &file : fs::directory_iterator(problemDir)) {
            if (file.is_regular_file() and file.path().filename() == ProblemFileName) {
                return file.path();
            }
        }

        throw std::runtime_error("no problem definition file found in '" + problemDir.string() + "'");
    }

    auto getProblems(const fs::path &problemDir) -> std::map<unsigned int, serialization::PartialProblem> {
        using namespace tempo::serialization;
        using namespace std::views;
        std::map<unsigned, PartialProblem> ret;
        const auto dir = problemDir / Serializer<>::SubProblemDir;
        if (not fs::is_directory(dir)) {
            throw std::runtime_error("no problems directory under " + problemDir.string());
        }

        const std::regex numberRegex("(\\d+)");
        for (const auto &file : fs::directory_iterator(problemDir / Serializer<>::SubProblemDir)) {
            if (file.is_regular_file() and
                file.path().filename().string().starts_with(Serializer<>::SubProblemBaseName)) {
                std::smatch res;
                const auto str = file.path().filename().string();
                if (std::regex_search(str, res, numberRegex) and res.size() == 2) {
                    auto id = std::stoi(*res.begin());
                    auto [_, inserted] = ret.emplace(id, deserializeFromFile<PartialProblem>(file));
                    if (not inserted) {
                        throw std::runtime_error("duplicate sub problem id " + std::to_string(id));
                    }
                }
            }
        }

        return ret;
    }

    auto loadDataPoint(const fs::path &mainDir, unsigned int id,
                       bool rootInstance) -> std::pair<DataPoint, DataPointStatus> {
        using namespace serialization;
        PartialProblem problem;
        if (rootInstance) {
            problem.associatedSolution = id;
        } else {
            auto problems = getProblems(mainDir);
            try {
                problem = problems.at(id);
            } catch (const std::out_of_range &) {
                return {{}, DataPointStatus::ProblemNotFound};
            }
        }

        Solution<int> solution;
        try {
            solution = tempo::getSolutions(mainDir).at(problem.associatedSolution);
        } catch (const std::out_of_range &) {
            return {{}, DataPointStatus::SolutionNotFound};
        }

        return {{std::move(problem), std::move(solution)}, DataPointStatus::Valid};
    }

    auto getSolutions(const fs::path &problemDir) -> std::map<unsigned int, serialization::Solution<int>> {
        using namespace tempo::serialization;
        std::map<unsigned, Solution<int>> ret;
        const auto dir = problemDir / Serializer<>::SolutionDir;
        if (not fs::is_directory(dir)) {
            throw std::runtime_error("no solutions directory under " + problemDir.string());
        }

        for (const auto &file : fs::directory_iterator(dir)) {
            if (file.is_regular_file() and
                file.path().filename().string().starts_with(Serializer<>::SolutionBaseName)) {
                auto sol = deserializeFromFile<Solution<int>>(file);
                auto id = sol.id;
                auto [_, inserted] = ret.emplace(id, std::move(sol));
                if (not inserted) {
                    throw std::runtime_error("duplicate sub problem id " + std::to_string(id));
                }
            }
        }

        return ret;
    }
}