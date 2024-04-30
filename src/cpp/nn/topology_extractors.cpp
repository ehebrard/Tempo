/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#include <Iterators.hpp>
#include <ranges>

#include "nn/tensor_utils.hpp"
#include "nn/topology_extractors.hpp"

namespace tempo::nn {

    MinimalTopologyBuilder::MinimalTopologyBuilder(std::size_t numEvents): edgeLookup(numEvents) {}

    MinimalTopologyBuilder::MinimalTopologyBuilder(const ProblemInstance &problem) : numTasks(problem.durations.size()),
                                                                                     numResources(
                                                                                             problem.resources.size()),
                                                                                     edgeLookup(numTasks * 2 + 2) {
        immutableTopologyFromDescription(problem);
    }

    auto MinimalTopologyBuilder::getTopology() -> Topology {
        return getTopology(makeSolverState(Matrix<int>{}));
    }

    void MinimalTopologyBuilder::resetToImmutableTopology() {
        edges.resize(immutableIdx);
        edgePairMask.resize(immutableIdx);
    }

    void MinimalTopologyBuilder::addEdge(const Edge &e, bool immutable, bool isResourceEdge, IndexType maskVal) {
        assert(not edgeLookup.contains(e));
        assert(not immutable or immutableIdx == edges.size());
        const auto idx = static_cast<IndexType>(edges.size());
        if (isResourceEdge) {
            edgeIdx.emplace_back(idx);
        }

        if (isResourceEdge or immutable) {
            edgeLookup(e) = idx;
            ++immutableIdx;
        }

        edges.emplace_back(e);
        edgePairMask.emplace_back(maskVal);
    }

    void MinimalTopologyBuilder::completeSubGraph(const Resource<int> &resourceSpec, IndexType resource) {
        const auto &tasks = resourceSpec;
        for (std::size_t i = 0; i < tasks.size(); ++i) {
            taskIdx.emplace_back(tasks[i]);
            resIdx.emplace_back(resource);
            resDemands.emplace_back(static_cast<DataType>(resourceSpec.getDemand(i)) /
                                    static_cast<DataType>(resourceSpec.capacity));
            for (std::size_t j = i + 1; j < tasks.size(); ++j) {
                Edge e(tasks[i], tasks[j]);
                Edge rev(tasks[j], tasks[i]);
                edgeRelResIdx.emplace_back(resource);
                edgeRelResIdx.emplace_back(resource);
                if (edgeLookup.contains(e)) {
                    edgeIdx.emplace_back(edgeLookup(e));
                    edgeIdx.emplace_back(edgeLookup(rev));
                    continue;
                }

                addEdge(e, true, true, pairMaskVal);
                addEdge(rev, true, true, pairMaskVal);
                ++pairMaskVal;
            }
        }
    }

    void MinimalTopologyBuilder::immutableTopologyFromDescription(const ProblemInstance &problem) {
        using namespace iterators;
        for (auto [r, resourceSpec]: enumerate(problem.resources, 0l)) {
            completeSubGraph(resourceSpec, r);
        }

        addPrecedenceEdges(problem.constraints, true);
        resourceDependencies = util::makeIndexTensor(taskIdx, resIdx);
        resourceDemands = torch::from_blob(resDemands.data(), {static_cast<long>(taskIdx.size()), 1},
                                           dataTensorOptions()).clone();
        edgeResourceRelations = util::makeIndexTensor(edgeIdx, edgeRelResIdx);
    }
}