/**
 * @author Tim Luchterhand
 * @date 25.07.23.
 */

#include <Iterators.hpp>

#include "nn/GNN.hpp"
#include "nn/tensor_utils.hpp"

namespace tempo::nn {
    GNN::GNN(const std::filesystem::path &modelLocation): model(torch::jit::load(modelLocation, Device)) {
        model.eval();
    }

    auto GNN::inference(const InputGraph &input) -> torch::IValue {
        return model.forward({input});
    }

    EdgeRegressor::EdgeRegressor(const std::filesystem::path &modelLocation) : GNN(modelLocation) {}

    auto EdgeRegressor::getHeatMap(const InputGraph &input) -> Matrix<DataType> {
        auto res = this->inference(input).toTensor().contiguous();
        return extractHeatMap(res, input.at(GraphKeys::EdgeIdx), input.at(GraphKeys::TaskFeatures).size(0));
    }

    auto EdgeRegressor::extractHeatMap(const torch::Tensor &edgeProbabilities, const torch::Tensor &edgeIdx,
                                       std::size_t numTasks) -> Matrix<DataType> {
        Matrix<DataType> ret(numTasks, numTasks, NoValue);
        const auto edgesFrom = util::getIndexSlice(edgeIdx, 0);
        const auto edgesTo = util::getIndexSlice(edgeIdx, 1);
        for (auto [idx, from, to] : iterators::zip_enumerate(edgesFrom, edgesTo, 0l)) {
            ret(from, to) = util::c_index<DataType>(edgeProbabilities, {idx});
        }

        return ret;
    }

    double EdgeRegressor::dstEdgeProbability(unsigned int taskFrom, unsigned int taskTo,
                                             const Matrix<DataType> &heatMap) {
        const auto edgeProb = heatMap(taskFrom, taskTo);
        const auto inverseProb = heatMap(taskTo, taskFrom);
        if (edgeProb == GNN::NoValue or inverseProb == GNN::NoValue) {
            std::stringstream ss;
            ss << "Invalid edge in GNN: task edge " << taskFrom << " -> " << taskTo << std::endl;
            throw std::runtime_error(ss.str());
        }

        // in accordance with https://arxiv.org/pdf/1605.02406.pdf based on DS evidence theory. For bernoulli
        // probabilities this is simply the edgeProbability
        return edgeProb + 0.5 * (1 - edgeProb - inverseProb);
    }
} // nn