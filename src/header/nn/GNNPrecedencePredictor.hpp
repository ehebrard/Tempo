/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#ifndef TEMPO_GNNPRECEDENCEPREDICTOR_HPP
#define TEMPO_GNNPRECEDENCEPREDICTOR_HPP

#include <vector>
#include <filesystem>
#include <ranges>
#include <Iterators.hpp>

#include "GNN.hpp"
#include "GraphBuilder.hpp"
#include "Literal.hpp"
#include "util/SchedulingProblemHelper.hpp"
#include "util/traits.hpp"
#include "torch_types.hpp"
#include "heat_map_utils.hpp"

namespace tempo {
    template<typename T>
    class Solver;
}

namespace tempo::nn {
    namespace fs = std::filesystem;

    template<concepts::scalar T, SchedulingResource R>
    class GNNPrecedencePredictor {
        EdgeRegressor gnn;
        GraphBuilder<T, R> graphBuilder;
        std::vector<Literal<T>> literals;
        std::vector<DataType> massesPos;
        std::vector<DataType> massesNeg;
    public:

        GNNPrecedencePredictor(const fs::path &modelLocation, const fs::path &featureExtractorConfigLocation,
                               SchedulingProblemHelper<T, R> problemInstance, std::vector<Literal<T>> literals) : gnn(
                modelLocation), graphBuilder(featureExtractorConfigLocation, std::move(problemInstance)), literals(
                std::move(literals)), massesPos(this->literals.size(), 0), massesNeg(massesPos) {
            for (auto l : literals) {
                if (not l.hasSemantic() or not l.isBoolean()) {
                    throw std::runtime_error("all literals need to be boolean and have a semantic");
                }
            }
        }

        void updateConfidence(const Solver<T> &solver) {
            const auto taskNetwork = graphBuilder.getProblem().getTaskDistances(solver);
            const auto graph = graphBuilder.getGraph(makeSolverState(std::move(taskNetwork), solver));
            const auto edgeHeatMap = gnn.getHeatMap(graph);
            const auto &mapping = graphBuilder.getProblem().getMapping();
            for (auto [lit, pos, neg] : iterators::zip(literals, massesPos, massesNeg)) {
                pos += probabilityMass(lit, edgeHeatMap, solver, mapping);
                neg += probabilityMass(~lit, edgeHeatMap, solver, mapping);
            }
        }

        auto getLiterals(double confidenceThreshold) const -> std::vector<Literal<T>> {
            std::vector<Literal<T>> ret;
            for (auto [lit, pos, neg] : iterators::zip(literals, massesPos, massesNeg)) {
                // = |m1 - m2| / (m1 + m2)
                auto certainty = std::abs(2 * pos / (pos + neg) - 1);
                if (certainty > confidenceThreshold) {
                    ret.emplace_back(pos > neg ? lit : ~lit);
                }
            }

            return ret;
        }

    };
}

#endif //TEMPO_GNNPRECEDENCEPREDICTOR_HPP
