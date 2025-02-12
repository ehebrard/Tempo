/**
* @author Tim Luchterhand
* @date 25.07.24
* @brief contains utility that allows to easily load and solve scheduling problems
*/

#include <ranges>
#include <Iterators.hpp>

#include "scheduling_helpers.hpp"
#include "util/SchedulingProblemHelper.hpp"


void loadBranch(tempo::Solver<int> &solver, const tempo::serialization::Branch &branch) {
    for (auto [id, val] : branch) {
        auto lit = solver.boolean.getLiteral(val, id);
        solver.set(lit);
    }

    solver.propagate();
}


EdgeMapper::EdgeMapper(const tempo::Options &options) {
    auto res = loadSchedulingProblem(options);
    solver = std::move(res.solver);
    problem = std::move(res.instance);
}

auto EdgeMapper::getTaskEdge(tempo::Literal<tempo::Time> lit) const -> std::pair<unsigned, unsigned> {
    if (not lit.hasSemantic()) {
        throw std::runtime_error("cannot get edge without semantic");
    }

    auto edge = solver->boolean.getEdge(lit);
    const auto &mapping = problem.getMapping();
    if (not mapping.contains(edge.from) or not mapping.contains(edge.to)) {
        throw std::runtime_error("edge does not correspond to task - task edge");
    }

    return {mapping(edge.from), mapping(edge.to)};
}

std::size_t EdgeMapper::numTasks() const noexcept {
    return problem.tasks().size();
}

auto EdgeMapper::getSolver() noexcept -> tempo::Solver<tempo::Time> & {
    return *solver;
}

auto EdgeMapper::getSolver() const noexcept -> const tempo::Solver<tempo::Time> & {
    return *solver;
}
