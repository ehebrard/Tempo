/**
 * @author Tim Luchterhand
 * @date 25.07.23.
 */

#include <Iterators.hpp>
#include <iostream>

#include "nn/GNN.hpp"
#include "nn/tensor_utils.hpp"

namespace tempo::nn {
    GNN::GNN(const std::filesystem::path &modelLocation): model(torch::jit::load(modelLocation, Device)) {
        model.eval();
        std::cout << "-- GNN CPU inference on " << torch::get_num_threads() << " threads" << std::endl;
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
} // nn
