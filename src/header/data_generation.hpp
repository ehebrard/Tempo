/**
* @author Tim Luchterhand
* @date 26.07.24
* @brief contains structures for training data generation
*/

#ifndef TEMPO_DATA_GENERATION_HPP
#define TEMPO_DATA_GENERATION_HPP

#include <filesystem>
#include <utility>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <vector>
#include <Iterators.hpp>
#include <map>
#include <nlohmann/json.hpp>

#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/serialization.hpp"
#include "util/TraceWatcher.hpp"

namespace tempo {
    namespace fs = std::filesystem;
    constexpr auto ProblemFileName = "problem_definition.txt";
    constexpr auto InfoFileName = "info.json";
    constexpr auto RootName = "root";
    constexpr auto GraphName = "graph";
    constexpr auto InputsName = "inputs";
    constexpr auto OutputsName = "outputs";
    constexpr auto LabelFileName = "task_network";
    constexpr auto LabelName = "label";
    constexpr auto GraphReferenceFile = "ref.json";

    /**
     * @brief Class that can be used to serialize solutions and partial problems to a files.
     * @details @copybrief
     * Automatically creates required folder structure and generates file names
     * @tparam T timing type
     */
    template<concepts::scalar T = int>
    class Serializer {
        fs::path targetDirectory;
        std::vector<serialization::Solution<T>> solutions;
        std::vector<serialization::PartialProblem> problems;
    public:
        static constexpr auto SolutionBaseName = "solution";
        static constexpr auto FileExtension = ".json";
        static constexpr auto SubProblemBaseName = "sub_problem";
        static constexpr auto SolutionDir = "solutions";
        static constexpr auto SubProblemDir = "sub_problems";

        /**
         * Ctor
         * @param targetDirectory directory where files should be saved
         */
        explicit Serializer(fs::path targetDirectory): targetDirectory(std::move(targetDirectory)) {
            std::filesystem::create_directories(this->targetDirectory / SolutionDir);
            std::filesystem::create_directory(this->targetDirectory / SubProblemDir);
        }

        /**
         * Adds a solution
         * @param objective objective value of the solution
         * @param decisions all decisions made by the solver
         */
        void addSolution(T objective, serialization::Branch decisions) {
            solutions.emplace_back(solutions.size(), objective, std::move(decisions));
        }

        /**
         * Adds a sub problem
         * @param decisions all decisions made by the solver
         */
        void addSubProblem(serialization::Branch decisions) {

#ifdef __DEBUG_BUILD__
            std::ranges::sort(decisions);
#endif
            problems.emplace_back(lastSolutionId(), std::move(decisions));
        }

        /**
         * Saves all solutions and sub problems to separate files under the target directory while respecting the folder
         * structure. Clears problem and solution buffers
         */
        void flush() {
            for (const auto &sol : solutions) {
                serializeToFile(sol, generateFileName(true, sol.id));
            }

            solutions.clear();
            for (auto [id, prob] : iterators::const_enumerate(problems)) {
                serializeToFile(prob, generateFileName(false, id));
            }

            problems.clear();
        }

        /**
         * Clears all registered solutions and problems
         */
        void clear() noexcept {
            solutions.clear();
            problems.clear();
        }

        /**
         * Gets the number of registered solutions
         * @return number of registered solutions
         */
        [[nodiscard]] std::size_t solutionCount() const noexcept {
            return solutions.size();
        }

        /**
         * Gets the number of registered sub problems
         * @return number of registered sub problems
         */
        [[nodiscard]] std::size_t problemCount() const noexcept {
            return problems.size();
        }

        /**
         * Access the solutions
         * @return const ref to solutions
         */
        [[nodiscard]] const auto &getSolutions() const noexcept {
            return solutions;
        }

        /**
         * Access the problems
         * @return const ref to problems
         */
        [[nodiscard]] const auto &getProblems() const noexcept {
            return problems;
        }

    private:
        [[nodiscard]] auto generateFileName(bool solution, unsigned id) const -> std::string {
            using namespace std::filesystem;
            std::stringstream ss;
            ss << (solution ? path(SolutionDir) / SolutionBaseName : path(SubProblemDir) / SubProblemBaseName).string()
               << id << FileExtension;
            return targetDirectory / ss.str();
        }

        [[nodiscard]] std::size_t lastSolutionId() const {
            if (solutions.size() == 0) {
                throw std::runtime_error("no solution has been found yet");
            }

            return solutions.size() - 1;
        }

    };

    /**
     * @brief Follows the solver's progress and generates solutions and sub problems at appropriate points.
     * @details @copybrief
     * @tparam T timing type
     */
    template<concepts::scalar T = int>
    class DataGenerator {
        Tracer tracer;
        Serializer<T> serializer;
        Interval<T> schedule;
        SubscriberHandle solutionHandler;
        SubscriberHandle deviationHandler;

    public:
        SubscribableEvent<const serialization::PartialProblem&, const serialization::Solution<T> &> DataPointCreated;
        ///< Triggers when a data point is found

        DataGenerator(DataGenerator &&) noexcept = default;
        DataGenerator &operator=(DataGenerator &&) noexcept = default;

        /**
         * Dtor. Saves all obtained solutions and problems
         */
        ~DataGenerator() {
            serializeData();
        }

        /**
         * Ctor
         * @param solver target solver
         * @param schedule interval representing the overall schedule
         * @param saveDirectory directory where to save the generated problems and solutions
         */
        DataGenerator(const Solver<T> &solver, Interval<T> schedule, fs::path saveDirectory) :
            tracer(solver), serializer(std::move(saveDirectory)), schedule(schedule),
            solutionHandler(solver.SolutionFound.subscribe_handled([this](const auto &solver) {
                this->handleSolution(solver);
            })),
            deviationHandler(tracer.DeviationOccurred.subscribe_handled(
                    [this](DeviationType type, const auto &conflicts, const auto &decisions) {
                switch (type) {
                    case DeviationType::Propagation:
                        handlePropagation(conflicts, decisions);
                        break;
                    case DeviationType::Fail:
                        handleConflict(decisions);
                        break;
                    case DeviationType::Decision:
                        break;
                    default:
                        throw std::runtime_error("enum out of bounds");
                }
            })) {}

        /**
         * Saves all generated solutions and problems under the save directory
         */
        void serializeData() {
            serializer.flush();
        }

        /**
         * Clears all found solutions and problems
         */
        void clear() noexcept {
            serializer.clear();
        }

        /**
         * Access the solutions
         * @return const ref to solutions
         */
        [[nodiscard]] std::size_t solutionCount() const noexcept {
            return serializer.solutionCount();
        }

        /**
         * Access the problems
         * @return const ref to problems
         */
        [[nodiscard]] std::size_t problemCount() const noexcept {
            return serializer.problemCount();
        }

    private:

        void handleConflict(serialization::Branch decisions) {
            serializer.addSubProblem(std::move(decisions));
            DataPointCreated.trigger(serializer.getProblems().back(), serializer.getSolutions().back());
        }

        void handlePropagation(const TraceWatcher::Conflicts &conflicts, serialization::Branch decisionsOnTrack) {
            for (auto [var, val] : conflicts) {
                decisionsOnTrack.emplace_back(var, val);
                serializer.addSubProblem(decisionsOnTrack);
                DataPointCreated.trigger(serializer.getProblems().back(), serializer.getSolutions().back());
                decisionsOnTrack.pop_back();
            }
        }

        void handleSolution(const Solver<T> &solver) {
            serialization::Branch decisions;
            const auto &traceWatcher = tracer.getWatcher();
            decisions.reserve(traceWatcher.getLastSolution().size());
            for (auto [var, val] : iterators::enumerate(traceWatcher.getLastSolution(), traceWatcher.getOffset())) {
                decisions.emplace_back(var, val);
            }

            serializer.addSolution(schedule.getEarliestEnd(solver), std::move(decisions));
        }

    };

    /**
     * Gets the path to the problem definition under a given problem directory
     * @param problemDir directory containing data points
     * @return path to instance file
     * @throws std::runtime_error if instance file could not be found under given path
     */
    auto getInstance(const fs::path &problemDir) -> fs::path;

    /**
     * Loads the problem info file
     * @param problemDir directory to the problem folder
     * @return deserialized problem info as json
     */
    auto getInfo(const fs::path &problemDir) -> nlohmann::json;

    /**
     * Loads all solutions under problem directory
     * @param problemDir path to directory with data points
     * @return map with deserialized solutions indexed by their ids
     * @throws std::runtime_error if solutions directory could not be found under given path
     */
    auto getSolutions(const fs::path &problemDir) -> std::map<unsigned int, serialization::Solution<int>>;

    /**
     * Loads all sub problems under problem directory
     * @param problemDir path to directory with data points
     * @return map with deserialized problems indexed by their ids
     * @throws std::runtime_error if sub_problems directory could not be found under given path
     */
    auto getProblems(const fs::path &problemDir) -> std::map<unsigned int, serialization::PartialProblem>;

    /**
     * @brief Datapoint load status
     * @details @copybrief
     */
    enum class DataPointStatus {
        Valid, ///< Data point was loaded correctly and is valid
        SolutionNotFound, ///< The associated solution to the data point was not found
        ProblemNotFound ///< Sub problem with the requested ID was not found
    };

    /**
     * @brief Represents a data point
     * @details @copybrief
     */
    struct DataPoint {
        serialization::PartialProblem problem; ///< partial problem representing the input
        serialization::Solution<int> solution; ///< expected locally optimal solution
    };

    /**
     * Loads a partial problem together with its associated solution
     * @param mainDir root directory of the problem
     * @param id number of the sub problem
     * @param rootInstance whether to load a root instance. In this case, the id is the id of the solution
     * @return DataPoint with problem and solution along with a status value.
     * Check the status to now if loading was successful
     */
    auto loadDataPoint(const fs::path &mainDir, unsigned id,
                       bool rootInstance) -> std::pair<DataPoint, DataPointStatus>;

}

#endif //TEMPO_DATA_GENERATION_HPP
