/**
 * @author Tim Luchterhand
 * @date 14.11.23.
 */

#include <stdexcept>
#include <sstream>

#include "nn/GNNEdgePolarityPredictor.hpp"
#include "Global.hpp"

namespace tempo::nn::heuristics {

    bool detail::choosePolarityFromHeatMap(unsigned taskFrom, unsigned taskTo, const Matrix<DataType> &heatMap) {
        const auto edgeProb = heatMap(taskFrom, taskTo);
        const auto inverseProb = heatMap(taskTo, taskFrom);
        if (edgeProb == GNN::NoValue or inverseProb == GNN::NoValue) {
            std::stringstream ss;
            ss << "Invalid edge in GNN: task edge " << taskFrom << " -> " << taskTo << std::endl;
            throw std::runtime_error(ss.str());
        }

        // in accordance with https://arxiv.org/pdf/1605.02406.pdf based on DS evidence theory. For bernoulli
        // probabilities this is simply the edgeProbability
        const auto pignisticProb = edgeProb + 0.5 * (1 - edgeProb - inverseProb);
        return pignisticProb > 0.5;

    }
}