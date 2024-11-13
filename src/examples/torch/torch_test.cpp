/**
 * @author Tim Luchterhand
 * @date 30.04.24.
 */
#include <iostream>

#include "nn/tensor_utils.hpp"

int main() {
    using namespace tempo::nn;
    auto tensor = util::makeIndexTensor({1, 2, 3}, {4, 5, 6});
    std::cout << "original tensor\n" << tensor << std::endl;
    auto slice = util::getIndexSlice(util::sliceAssign(tensor, {1, util::SliceHere}, {-1l, -2l, -3l}), 1);
    std::cout << "modified lower slice" << std::endl;
    for (auto v : slice) {
        std::cout << v << " ";
    }

    std::cout << std::endl;
    std::cout << "modified tensor\n" << tensor << std::endl;
    std::cout << "edges view" << std::endl;
    for (auto &&[f, t] : util::getEdgeView(tensor)) {
        std::cout << f << " -> " << t << "\n";
    }

    std::cout << std::endl;
    return 0;
}