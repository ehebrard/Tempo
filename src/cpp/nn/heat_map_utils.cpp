/**
* @author Tim Luchterhand
* @date 11.09.24
* @brief
*/

#include "nn/heat_map_utils.hpp"

namespace tempo::nn {

    DataType probabilityMass(unsigned int taskFrom, unsigned int taskTo, const Matrix<DataType> &heatMap) {
        auto prob = heatMap(taskFrom, taskTo);
        if (prob == GNN::NoValue) {
            std::stringstream ss;
            ss << "Invalid edge in GNN: task edge " << taskFrom << " -> " << taskTo << std::endl;
            throw std::runtime_error(ss.str());
        }

        return prob;
    }
}