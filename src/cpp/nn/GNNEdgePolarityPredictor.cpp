/**
 * @author Tim Luchterhand
 * @date 14.11.23.
 */

#include <stdexcept>
#include <sstream>

#include "nn/GNNEdgePolarityPredictor.hpp"
#include "Global.hpp"

namespace tempo::nn::heuristics {

    GNNEdgePolarityPredictor::GNNEdgePolarityPredictor(const std::filesystem::path &modelLocation,
                                                       const std::filesystem::path &featureExtractorConfigLocation,
                                                       const ProblemInstance &problemInstance) :
            gnn(modelLocation),
            graphBuilder(featureExtractorConfigLocation, problemInstance) {}

    bool GNNEdgePolarityPredictor::choosePolarityFromHeatMap(event from, event to, const Matrix<DataType> &heatMap) {
        from = TASK(from);
        to = TASK(to);
        const auto edgeProb = heatMap(from, to);
        const auto inverseProb = heatMap(to, from);
        if (edgeProb == GNN::NoValue or inverseProb == GNN::NoValue) {
            std::stringstream ss;
            ss << "Invalid edge in GNN: task edge " << from << " -> " << to << std::endl;
            throw std::runtime_error(ss.str());
        }

        // in accordance with https://arxiv.org/pdf/1605.02406.pdf based on DS evidence theory. For bernoulli
        // probabilities this is simply the edgeProbability
        const auto pignisticProb = edgeProb + 0.5 * (1 - edgeProb - inverseProb);
        return pignisticProb > 0.5;

    }
}