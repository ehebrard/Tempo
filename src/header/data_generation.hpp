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

#include "util/traits.hpp"
#include "util/SubscribableEvent.hpp"
#include "util/serialization.hpp"
#include "util/TraceWatcher.hpp"

namespace tempo {
    namespace fs = std::filesystem;

    /**
     * @brief Class that can be used to serialize solutions and partial problems to a files.
     * @details @copybrief
     * Automatically creates required folder structure and generates file names
     * @tparam T timing type
     */
    template<concepts::scalar T = int>
    class Serializer {
        fs::path targetDirectory;
        unsigned solutionId = 0;
        unsigned problemId = 0;
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
            solutions.emplace_back(solutionId++, objective, std::move(decisions));
        }

        /**
         * Adds a sub problem
         * @param decisions all decisions made by the solver
         */
        void addSubProblem(serialization::Branch decisions) {
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
            for (const auto &prob : problems) {
                serializeToFile(prob, generateFileName(false, problemId++));
            }

            problems.clear();
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
            if (solutionId == 0) {
                throw std::runtime_error("no solution has been found yet");
            }

            return solutionId - 1;
        }

    };

    /**
     * @brief Follows the solver's progress and generates solutions and sub problems at appropriate points.
     * @details @copybrief
     * @tparam T timing type
     */
    template<concepts::scalar T = int>
    class Tracer {
        class Aligned {
            const Solver<T> &s;
        public:
            explicit constexpr Aligned(const Solver<T> &solver) noexcept: s(solver) {}
            constexpr bool operator()(var_t var, bool truthVal) const noexcept {
                return s.boolean.isUndefined(var) or s.boolean.isTrue(var) == truthVal;
            }
        };

        TraceWatcher traceWatcher;
        Serializer<T> serializer;
        Interval<T> schedule;
        SubscriberHandle solutionHandler;
        SubscriberHandle decisionHandler;
        SubscriberHandle conflictHandler;
        SubscriberHandle backtrackHandler;
        SubscriberHandle propagationHandler;

    public:
        Tracer(Tracer &&) noexcept = default;
        Tracer &operator=(Tracer &&) noexcept = default;

        /**
         * Dtor. Saves all obtained solutions and problems
         */
        ~Tracer() {
            serializeData();
        }

        /**
         * Ctor
         * @param solver target solver
         * @param schedule interval representing the overall schedule
         * @param saveDirectory directory where to save the generated problems and solutions
         */
        Tracer(const Solver<T> &solver, Interval<T> schedule, fs::path saveDirectory) :
            traceWatcher(solver.boolean_search_vars), serializer(std::move(saveDirectory)), schedule(schedule),
            solutionHandler(solver.SolutionFound.subscribe_handled([this](const auto &solver) {
                this->handleSolution(solver);
            })),
            decisionHandler(solver.ChoicePoint.subscribe_handled([this](auto lit) {
                this->traceWatcher.step(lit);
            })),
            conflictHandler(solver.ConflictEncountered.subscribe_handled([this, &solver](const auto &) {
                this->handleConflict(solver);
            })),
            backtrackHandler(solver.BackTrackCompleted.subscribe_handled([this, &solver]() {
                this->traceWatcher.updateOnTrack(Aligned(solver));
            })),
            propagationHandler(solver.PropagationCompleted.subscribe_handled([this](const auto &solver) {
                this->handlePropagation(solver);
            })) {}

        /**
         * Saves all generated solutions and problems under the save directory
         */
        void serializeData() {
            serializer.flush();
        }

    private:
        static auto getDecisions(const Solver<T> &solver) -> serialization::Branch {
            serialization::Branch branch;
            std::ranges::subrange currentDecisions(solver.getBranch().bbegin(), solver.getBranch().bend());
            branch.reserve(currentDecisions.size() + 1);
            for (var_t var : currentDecisions) {
                assert(not solver.boolean.isUndefined(var));
                branch.emplace_back(var, solver.boolean.isTrue(var));
            }

            return branch;
        }

        void handleConflict(const Solver<T> &solver) {
            if (not traceWatcher.isOnTrack()) {
                return;
            }

            serializer.addSubProblem(getDecisions(solver));
            traceWatcher.setOnTrack(false);
        }

        void handlePropagation(const Solver<T> &solver) {
            if (not traceWatcher.isOnTrack()) {
                return;
            }

            const auto conflicting = traceWatcher.updateOnTrack(Aligned(solver));
            auto currentDecisions = getDecisions(solver);
            for (auto [var, val] : conflicting) {
                currentDecisions.emplace_back(var, val);
                serializer.addSubProblem(currentDecisions);
                currentDecisions.pop_back();
            }
        }

        void handleSolution(const Solver<T> &solver) {
            traceWatcher.registerSolution([&solver](auto var) {
                assert(not solver.boolean.isUndefined(var));
                return solver.boolean.isTrue(var);
            });

            serialization::Branch decisions;
            decisions.reserve(traceWatcher.getLastSolution().size());
            for (auto [var, val] : iterators::enumerate(traceWatcher.getLastSolution(), traceWatcher.getOffset())) {
                decisions.emplace_back(var, val);
            }

            serializer.addSolution(schedule.getEarliestEnd(solver), std::move(decisions));
        }

    };

}

#endif //TEMPO_DATA_GENERATION_HPP
