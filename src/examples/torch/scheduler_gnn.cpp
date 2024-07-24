/**
* @author Tim Luchterhand
* @date 10.07.24
* @brief
*/

#include "Solver.hpp"
#include "heuristics/GNNValueHeuristics.hpp"
#include "heuristics/heuristic_factories.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/osp.hpp"


int main(int argc, char **argv) {
    using namespace tempo;
    using namespace heuristics;
    Parser p = getBaseParser();
    std::string gnnLocation;
    std::string featureExtractorConf;
    p.getCmdLine().add<TCLAP::ValueArg<std::string>>(gnnLocation, "", "gnn-loc", "Location of the trained GNN", true,
                                                     "", "string");
    p.getCmdLine().add<TCLAP::ValueArg<std::string>>(featureExtractorConf, "", "feat-config",
                                                     "Location of the feature extractor config", true, "", "string");
    p.parse(argc, argv);
    const auto &opt = p.getOptions();
    Solver<> S(opt);
    auto schedule{S.newInterval(0, Constant::Infinity<int>, 0, 0, 0, Constant::Infinity<int>)};
    std::vector<DisjunctiveResource<>> resources;
    std::vector<Interval<>> tasks;

    if (opt.input_format == "osp") {
        osp::parse(opt.instance_file, S, schedule, tasks, resources);
    } else {
        std::cerr << "problem type " << opt.input_format << " is not (yet) supported" << std::endl;
        std::exit(1);
    }

    SchedulingProblemHelper problem(std::move(tasks), std::move(resources), {}, schedule);

    for (const auto &consumingTasks: problem.resources()) {
        auto constraint = NoOverlap(schedule, consumingTasks);
        S.post(constraint);
    }

    auto trivialUb = 0;
    for (const auto &t : problem.tasks()) {
        trivialUb += t.minDuration(S);
    }

    trivialUb = std::min(trivialUb, opt.ub);
    S.set(schedule.end.before(trivialUb));

    GNNFullGuidance valBranching(opt.polarity_epsilon, gnnLocation, featureExtractorConf, std::move(problem));
    auto varBranching = make_variable_heuristic(S);
    S.setBranchingHeuristic(make_compound_heuristic(std::move(varBranching), std::move(valBranching)));
    S.minimize(schedule.end);
    return 0;
}