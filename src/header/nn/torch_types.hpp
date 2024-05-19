//
// Created by tim on 20.07.23.
//

#ifndef SCHEDCL_TORCH_TYPES_HPP
#define SCHEDCL_TORCH_TYPES_HPP

#include <c10/core/TensorOptions.h>
#include <ATen/core/Dict.h>
#include <torch/types.h>
#include <vector>


namespace tempo::nn {
#ifdef TORCH_USE_GPU
    constexpr auto Device = torch::kCUDA;
#else
    constexpr auto Device = at::kCPU;
#endif
    using IndexType = long;
    using DataType = float;
    using Edge = std::pair<IndexType, IndexType>;
    using EdgeVector = std::vector<Edge>;
    using InputGraph = c10::Dict<std::string, torch::Tensor>;

    /**
     * Tensor options for index tensors
     * @return
     */
    inline auto indexTensorOptions() noexcept -> torch::TensorOptions {
        return torch::dtype<IndexType>().memory_format(torch::MemoryFormat::Contiguous).device(Device);
    }

    /**
     * Tensor options for data tensors
     * @return
     */
    inline auto dataTensorOptions() noexcept -> torch::TensorOptions {
        return torch::dtype<DataType>().memory_format(torch::MemoryFormat::Contiguous).device(Device);
    }

    /**
     * Struct representing the graph topology corresponding to the solver's state
     */
    struct Topology {
        unsigned numTasks; ///< number of tasks in the problem
        unsigned numResources; ///< number of resources in the problem
        torch::Tensor edgeIndices; ///< all edges in COO format: [2 x num edges], [0:] src edges, [1:] destination edges
        torch::Tensor edgePairMask; ///< mask marking edge pairs: [num edges]. corresponding edges have the same mask value
        torch::Tensor edgeResourceRelations; ///< edges belonging to certain resources in COO format:
                                             ///< [2 x num edge resource relations], [0:] edge idx, [1:] res idx
        torch::Tensor resourceDependencies; ///< dependencies between tasks and resources in COO format:
                                            ///< [2 x num resource deps], [0:] task idx, [1:] resource idx
        torch::Tensor resourceDemands; ///< normalized resource dependencies: [num resource deps, 1]
    };

    /**
     * Defines dictionary keys for different graph properties and features
     */
    struct GraphKeys {
        static constexpr auto TaskFeatures = "task_features"; ///< name of the task feature tensor
        static constexpr auto EdgeFeatures = "edge_features"; ///< name of the edge feature tensor
        static constexpr auto ResourceFeatures = "resource_features"; ///< name of the resource feature tensor
        static constexpr auto ResourceConsumptions = "resource_consumptions"; ///< name of the resource consumptions tensor
        static constexpr auto EdgeIdx = "edge_idx"; ///< name of the edge index tensor
        static constexpr auto ResourceDependencies = "resource_dependencies"; ///< name of the resource dependencies index tensor
        static constexpr auto EdgeResourceRelations = "edge_resource_relations"; ///< name of the edge-resource-relation index tensor
        static constexpr auto EdgePairMask = "edge_pair_mask"; ///< name of the edge-pair mask
        static constexpr std::array AllKeys{TaskFeatures, EdgeFeatures, ResourceFeatures, ResourceConsumptions, EdgeIdx,
                                            ResourceDependencies, EdgeResourceRelations, EdgePairMask};
    };


}

#endif //SCHEDCL_TORCH_TYPES_HPP
