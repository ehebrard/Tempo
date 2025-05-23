/**
* @author Tim Luchterhand
* @date 06.08.24
* @brief Calculates the importance of binary decisions for all search variables in a given (sub) problem
*/

#include <string>
#include <limits>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include <optional>
#include <future>
#include <tuple>
#include <chrono>

#include "helpers/shell.hpp"
#include "helpers/cli.hpp"
#include "helpers/git_sha.hpp"
#include "helpers/scheduling_helpers.hpp"
#include "data_generation.hpp"
#include "util/KillHandler.hpp"
#include "helpers/VarImportanceRunner.hpp"
#include "heuristics/LNS/PolicyDecay.hpp"
#include "heuristics/LNS/relaxation_policy_factories.hpp"
#include "util/ThreadPool.hpp"

#define JSONIFY(JSON, DATA) JSON[#DATA] = DATA;

namespace fs = std::filesystem;

class RunResult {
public:
    using Run = std::optional<std::pair<Result, Result>>;
private:
    tempo::Literal<Time> lit;
    std::future<Run> aResult;
    Run sResult;
    bool synchronous;
public:
    RunResult(tempo::Literal<Time> lit, Run result) noexcept: lit(lit), sResult(std::move(result)), synchronous(true) {}

    RunResult(tempo::Literal<Time> lit, std::future<Run> result) noexcept: lit(lit), aResult(std::move(result)),
                                                                           synchronous(false) {}

    template<typename Dur>
    bool wait(Dur dur) {
        return sResult.has_value() or aResult.wait_for(dur) == std::future_status::ready;
    }

    auto get() {
        if (synchronous) {
            return std::make_pair(lit, sResult);
        }

        if (aResult.valid()) {
            return std::make_pair(lit, aResult.get());
        }

        throw std::runtime_error("No run result");
    }
};

int main(int argc, char **argv) {
    using namespace std::chrono_literals;
    using namespace tempo;
    using namespace heuristics;
    int subId = -1;
    bool root;
    bool useLns = false;
    lns::RelaxationPolicyParams policyParams{
        .decayConfig = lns::PolicyDecayConfig(), .numScheduleSlices = 4, .allTaskEdges = false
    };
    auto policyType = lns::RelaxPolicy::RandomTasks;
    double sporadicIncrement = 0;
    unsigned mt = 0;
    auto opt = cli::parseOptions(argc, argv,
                                 cli::ArgSpec("sub-number", "id of the sub problem", false, subId),
                                 cli::SwitchSpec("root", "calculate edge importance for a root instance instead", root,
                                                 false),
                                 cli::SwitchSpec("lns", "whether to use LNS", useLns, false),
                                 cli::ArgSpec("fix-decay", "relaxation ratio decay",
                                              false, policyParams.decayConfig.decay, 0.5),
                                 cli::ArgSpec("fix-ratio", "initial relaxation ratio",
                                              false, policyParams.decayConfig.fixRatio, 0.5),
                                 cli::ArgSpec("decay-min-fail",
                                              "lower bound solver failure rate for ratio decay config",
                                              false, policyParams.decayConfig.minFailRatio),
                                 cli::ArgSpec("decay-max-fail",
                                              "upper bound solver failure rate for ratio decay config",
                                              false, policyParams.decayConfig.maxFailRatio),
                                 cli::SwitchSpec("decay-on-success", "whether to decrease fix rate even on success",
                                                 policyParams.decayConfig.decreaseOnSuccess, false),
                                 cli::ArgSpec("retry-limit", "number of fails before decreasing relaxation ratio",
                                              false, policyParams.decayConfig.retryLimit),
                                 cli::ArgSpec("decay-mode", "relaxation ratio decay mode on failure", false,
                                              policyParams.decayConfig.decayMode),
                                 cli::ArgSpec("relax-slices", "number of schedule slices",
                                              false, policyParams.numScheduleSlices, 4),
                                 cli::SwitchSpec("fix-all-task-edges",
                                                 "whether to fix all task edges or only those between fixed tasks",
                                                 policyParams.allTaskEdges, false),
                                 cli::ArgSpec("lns-policy", "lns relaxation policy", false, policyType),
                                 cli::ArgSpec("sporadic-increment", "sporadic root search probability increment", false,
                                              sporadicIncrement),
                                 cli::ArgSpec("mt", "multi-threading number of threads", false, mt)
    );
    const auto mainDir = opt.instance_file;
    opt.instance_file = getInstance(mainDir);
    DataPoint dataPoint;
    if (root) {
        auto solutions = getSolutions(mainDir);
        if (solutions.empty()) {
            std::cerr << "No solutions found" << std::endl;
            std::exit(1);
        }

        dataPoint.solution = std::move(solutions.rbegin()->second);
    } else {
        DataPointStatus status;
        std::tie(dataPoint, status) = loadDataPoint(mainDir, subId, false);
        if (status == DataPointStatus::ProblemNotFound) {
            std::cerr << "sub problem with id " << subId << " could not be found" << std::endl;
            std::exit(1);
        } else if (status == DataPointStatus::SolutionNotFound) {
            std::cerr << "associated solution" << subId << " could not be found" << std::endl;
            std::exit(1);
        }
    }

    VarImportanceRunner runner(dataPoint.problem, opt, dataPoint.solution);
    if (runner.isInconsistent()) {
        std::cout << "snapshot is inconsistent" << std::endl;
        std::exit(2);
    }

    unsigned numIndifferent = 0;
    unsigned numDecisions = 0;
    unsigned numDecided = 0;
    unsigned numEdges = 0;
    const auto optimum = dataPoint.solution.objective;
    std::unordered_map<var_t, double> varImportance;
    const EdgeMapper mapper(opt);
    std::vector<std::string> taskEdgeImportance;
    std::optional<lns::RelaxationPolicy<Time, ResourceConstraint>> lnsBasePolicy;
    if (useLns) {
        auto problem = loadSchedulingProblem(opt);
        lnsBasePolicy = lns::make_relaxation_policy(policyType, problem.instance.tasks(), problem.constraints,
                                                    policyParams, opt.verbosity);
    }

    auto run = [useLns, sporadicIncrement, &lnsBasePolicy](auto &s, auto &obj) {
        if (useLns) {
            auto policyCopy = *lnsBasePolicy;
            auto policy = lns::make_sporadic_root_search(sporadicIncrement, policyCopy);
            s.largeNeighborhoodSearch(obj, policy);
        } else {
            s.optimize(obj);
        }
    };

    std::vector<RunResult> runs;
    ThreadPool threadPool(mt);
    std::mutex mutex;
    unsigned iteration = 0;
    unsigned progress = 0;
    for (auto lit: runner.getLiterals()) {
        if (varImportance.contains(lit.variable())) {
            //This should never happen
            std::cout << "duplicate literal in search variables" << std::endl;
            std::exit(2);
        }

        if (KillHandler::instance().signalReceived()) {
            break;
        }

        auto work = [&runner, lit, &run, &mutex, &iteration, &progress]() -> RunResult::Run {
            using enum SchedulerState;
            auto posResult = runner.run(lit, run);

            if (posResult.state == AlreadyDecided) {
                return {};
            }

            auto negResult = runner.run(~lit, run);
            if (negResult.state == AlreadyDecided) {
                return {};
            }

            std::lock_guard lock(mutex);
            const auto percentage = static_cast<unsigned>(
                static_cast<double>(++iteration) / runner.getLiterals().size() * 100);
            if (percentage > progress) {
                progress = percentage;
                std::cout << "-- " << progress << "% done" << std::endl;
            }
            return std::make_pair(posResult, negResult);
        };

        if (mt > 1) {
            runs.emplace_back(lit, threadPool.submit(work));
        } else {
            runs.emplace_back(lit, work());
        }
    }

    auto it = runs.begin();
    while (not runs.empty()) {
        using enum SchedulerState;
        auto &runResult = *it;
        if (not runResult.wait(1s)) {
            ++it;
            if (it == runs.end()) {
                it = runs.begin();
            }

            continue;
        }

        ++numEdges;
        auto [lit, res] = runResult.get();
        if (not res.has_value()) {
            ++numDecided;
            continue;
        }

        ++numDecisions;
        auto [posResult, negResult] = *res;
        if (posResult.result.value_or(0) != optimum and negResult.result.value_or(0) != optimum) {
            std::cout << "ERROR: " << lit << " (variable " << lit.variable() << ")"
                      << " is a choice point but the optimal solution cannot be found from this state" << std::endl;
            std::exit(2);
        }

        if (posResult.state == Unsat and negResult.state == Unsat) {
            std::cout << "ERROR: " << lit << " (variable " << lit.variable() << ")"
                      << "cannot be set! Inconsistent datapoint" << std::endl;
            std::exit(2);
        }

        constexpr auto Inf = std::numeric_limits<double>::infinity();
        const auto resPos = static_cast<std::optional<double>>(posResult.result).value_or(Inf);
        const auto resNeg = static_cast<std::optional<double>>(negResult.result).value_or(Inf);
        const auto larger = std::max(resPos, resNeg);
        const auto smaller = std::min(resPos, resNeg);
        const auto quality = 1 - static_cast<double>(smaller) / larger;
        numIndifferent += quality == 0;
        varImportance[lit.variable()] = quality;
        auto [f, t] = mapper.getTaskEdge(lit);
        taskEdgeImportance.emplace_back(
                std::to_string(f) + " <-> " + std::to_string(t) + "=" + std::to_string(quality));
        if (it == runs.end() - 1) {
            runs.pop_back();
            it = runs.begin();
        } else {
            std::swap(*it, runs.back());
            runs.pop_back();
        }
    }

    auto unimportant = static_cast<double>(numIndifferent) / numDecisions;
    auto decided = static_cast<double>(numDecided) / numEdges;
    if (std::isnan(unimportant)) { unimportant = 0; }
    if (std::isnan(decided)) { decided = 1; }
    nlohmann::json json;
    JSONIFY(json, varImportance);
    JSONIFY(json, taskEdgeImportance);
    JSONIFY(json, unimportant);
    JSONIFY(json, decided);
    json["problem"] = fs::path(mainDir).filename();
    json["searchEffort"] = runner.averageSearchEffort();
    json["avgDecisions"] = runner.averageNumberOfDecisions();
    json["date"] = shell::getTimeStamp();
    json["commit"] = GitSha;
    std::cout << json.dump(__JSON_INDENT__) << std::endl;
    return 0;
}
