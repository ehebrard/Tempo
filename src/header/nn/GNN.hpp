/**
 * @author Tim Luchterhand
 * @data 25.07.23.
 */

#ifndef TEMPO_GNN_HPP
#define TEMPO_GNN_HPP
#include <filesystem>
#include <torch/script.h>
#include <torchscatter/scatter.h>
#include <torchsparse/sparse.h>

#include "torch_types.hpp"
#include "util/Matrix.hpp"

namespace tempo::nn {

    /**
     * Graph neural network torch wrapper
     */
    class GNN {
    public:
        static constexpr DataType NoValue = -1;
        /**
         * Ctor
         * @param modelLocation path to the torch::jit model checkpoint
         */
        explicit GNN(const std::filesystem::path &modelLocation);

        /**
         * Runs the inference
         * @param input Input graph containing all features and the graph topology
         * @return output of the model
         */
        auto inference(const InputGraph &input) -> torch::IValue;

    protected:
        torch::jit::Module model;
    };

    /**
     * GNN that computes probabilities for each edge (an edge heat map)
     */
    class EdgeRegressor : public GNN {
    public:
        /**
         * CTor
         * @param modelLocation path to the torch::jit model checkpoint
         */
        explicit EdgeRegressor(const std::filesystem::path &modelLocation);

        /**
         * Runs the inference and constructs a matrix containing the edge probabilities
         * @param input input graph containing all features and the topology
         * @return task network matrix where (t, u) = probability that task t is scheduled after task u
         */
        auto getHeatMap(const InputGraph &input) -> Matrix<DataType>;

        /**
         * Computes the bayesian probability of an edge from taskFrom to taskTo using Dempster-Shafer theory
         * @param taskFrom source node of the edge
         * @param taskTo destination node of the edge
         * @param heatMap heat map with arbitrary probability masses
         * @return probability value in [0, 1]
         */
        static double dstEdgeProbability(unsigned taskFrom, unsigned taskTo, const Matrix<DataType> &heatMap);
    protected:
        static auto extractHeatMap(const torch::Tensor &edgeProbabilities, const torch::Tensor &edgeIdx,
                            std::size_t numTasks) -> Matrix<DataType>;
        using GNN::inference;
    };

} // nn

#endif //TEMPO_GNN_HPP
