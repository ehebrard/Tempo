/**
* @author Tim Luchterhand
* @date 30.07.24
* @brief
*/

#include <regex>
#include <functional>
#include <ranges>

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
        using namespace std::views;
        std::vector<PartialProblem> ret;
        const auto dir = problemDir / Serializer<>::SubProblemDir;
        if (not fs::is_directory(dir)) {
            throw std::runtime_error("no problems directory under " + problemDir.string());
        }

        std::vector<unsigned> ids;
        const std::regex numberRegex("(\\d+)");
        for (const auto &file : fs::directory_iterator(problemDir / Serializer<>::SubProblemDir)) {
            if (file.is_regular_file() and
                file.path().filename().string().starts_with(Serializer<>::SubProblemBaseName)) {
                std::smatch res;
                const auto str = file.path().filename().string();
                if (std::regex_search(str, res, numberRegex)) {
                    if (res.size() == 2) {
                        ret.emplace_back(deserializeFromFile<PartialProblem>(file));
                        ids.emplace_back(std::stoi(*res.begin()));
                    }
                }
            }
        }

        std::vector<std::pair<std::reference_wrapper<PartialProblem>, unsigned>> sorter;
        for (auto [p, id] : iterators::zip(ret, ids)) {
            sorter.emplace_back(p, id);
        }

        std::ranges::sort(sorter, [](const auto &a, const auto &b) { return std::get<1>(a) < std::get<1>(b); });
        auto elems = sorter | elements<0> | transform([](auto &pref) -> auto & { return pref.get(); });
        ret = {std::move_iterator(elems.begin()), std::move_iterator(elems.end())};
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

    auto getSolutions(const fs::path &problemDir) -> std::vector<serialization::Solution<int>> {
        using namespace tempo::serialization;
        std::vector<Solution<int>> ret;
        const auto dir = problemDir / Serializer<>::SolutionDir;
        if (not fs::is_directory(dir)) {
            throw std::runtime_error("no solutions directory under " + problemDir.string());
        }

        for (const auto &file : fs::directory_iterator(dir)) {
            if (file.is_regular_file() and
                file.path().filename().string().starts_with(Serializer<>::SolutionBaseName)) {
                ret.emplace_back(deserializeFromFile<Solution<int>>(file));
            }
        }

        std::ranges::sort(ret, [](const auto &a, const auto &b) { return a.id < b.id; });
        return ret;
    }
}