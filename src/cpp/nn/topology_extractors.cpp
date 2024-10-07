/**
 * @author Tim Luchterhand
 * @date 17.11.23.
 */

#include <Iterators.hpp>
#include <ranges>

#include "nn/tensor_utils.hpp"
#include "nn/topology_extractors.hpp"

namespace tempo::nn {

    auto MinimalTopologyBuilder::getTopology() const -> const Topology& {
        return cache;
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

    std::size_t impl::EdgeLookup::getNumEdges() const {
        std::size_t ret = 0;
        this->for_each([&ret](auto val) { ret += (val != EdgeLookup::NoValue); });
        return ret;
    }
}