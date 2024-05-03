/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#include <Iterators.hpp>
#include <ranges>

#include "nn/tensor_utils.hpp"
#include "nn/topology_extractors.hpp"

namespace tempo::nn {

    MinimalTopologyBuilder::MinimalTopologyBuilder(const ProblemInstance &problem) {
        using namespace iterators;
        cache.numTasks = problem.durations.size();
        cache.numResources = problem.resources.size();
        impl::TopologyData topologyData{.edgeLookup = impl::EdgeLookup(cache.numTasks * 2 + 2)};
        for (auto [r, resourceSpec]: enumerate(problem.resources, 0l)) {
            completeSubGraph(resourceSpec, r, topologyData);
        }

        addPrecedenceEdges(problem.constraints, topologyData);
        cache.edgeIndices = util::makeIndexTensor(topologyData.edges);
        cache.edgePairMask = torch::from_blob(topologyData.edgePairMask.data(),
                                              static_cast<long>(topologyData.edgePairMask.size()),
                                              indexTensorOptions()).clone();
        cache.edgeResourceRelations = util::makeIndexTensor(topologyData.edgeIdx, topologyData.edgeRelResIdx);
        cache.resourceDependencies = util::makeIndexTensor(topologyData.taskIdx, topologyData.resIdx);
        cache.resourceDemands = torch::from_blob(topologyData.resDemands.data(),
                                           {static_cast<long>(topologyData.taskIdx.size()), 1},
                                           dataTensorOptions()).clone();
    }

    auto MinimalTopologyBuilder::getTopology() const -> const Topology& {
        return getTopology(makeSolverState(Matrix<int>{}));
    }

    void MinimalTopologyBuilder::addEdge(const Edge &e, bool isResourceEdge, IndexType maskVal, impl::TopologyData &topologyData) {
        assert(not topologyData.edgeLookup.contains(e));
        const auto idx = static_cast<IndexType>(topologyData.edges.size());
        if (isResourceEdge) {
            topologyData.edgeIdx.emplace_back(idx);
        }

        topologyData.edgeLookup(e) = idx;
        topologyData.edges.emplace_back(e);
        topologyData.edgePairMask.emplace_back(maskVal);
    }

    void MinimalTopologyBuilder::completeSubGraph(const Resource<int> &resourceSpec, IndexType resource,
                                                  impl::TopologyData &topologyData) {
        const auto &tasks = resourceSpec;
        for (std::size_t i = 0; i < tasks.size(); ++i) {
            topologyData.taskIdx.emplace_back(tasks[i]);
            topologyData.resIdx.emplace_back(resource);
            topologyData.resDemands.emplace_back(static_cast<DataType>(resourceSpec.getDemand(i)) /
                                                 static_cast<DataType>(resourceSpec.capacity));
            for (std::size_t j = i + 1; j < tasks.size(); ++j) {
                Edge e(tasks[i], tasks[j]);
                Edge rev(tasks[j], tasks[i]);
                topologyData.edgeRelResIdx.emplace_back(resource);
                topologyData.edgeRelResIdx.emplace_back(resource);
                if (topologyData.edgeLookup.contains(e)) {
                    topologyData.edgeIdx.emplace_back(topologyData.edgeLookup(e));
                    topologyData.edgeIdx.emplace_back(topologyData.edgeLookup(rev));
                    continue;
                }

                addEdge(e, true, topologyData.pairMaskVal, topologyData);
                addEdge(rev, true, topologyData.pairMaskVal, topologyData);
                ++topologyData.pairMaskVal;
            }
        }
    }

    std::size_t impl::EdgeLookup::getNumEdges() const {
        std::size_t ret = 0;
        this->for_each([&ret](auto val) { ret += (val != EdgeLookup::NoValue); });
        return ret;
    }
}