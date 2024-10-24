/**
* @author Tim Luchterhand
* @date 06.08.24
* @brief simple program that loads two input graphs and checks if they're equivalent
*/

#include <iostream>

#include "nn/tensor_utils.hpp"

int main(int argc, char **argv) {
    using namespace tempo::nn::util;
    if (argc != 3) {
        std::cerr << "specify the paths to the two graphs" << std::endl;
        std::exit(1);
    }

    auto a = loadGraph(argv[1]);
    auto b = loadGraph(argv[2]);
    if (not compareGraphs(a, b)) {
        std::cout << "graphs not equivalent" << std::endl;
        std::exit(-1);
    }

    std::cout << "graphs equivalent" << std::endl;
    return 0;
}