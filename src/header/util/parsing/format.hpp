#ifndef __FORMAT_HH
#define __FORMAT_HH

#include <vector>

#include "Global.hpp"


struct ProblemInstance {
    int lowerBound;
    int optimalSolution;
    std::vector<int> durations;
    std::vector<std::tuple<tempo::event, tempo::event, int>> constraints;
    std::vector<std::vector<tempo::task>> resources;
    //  std::vector<std::vector<int>> jobs;
};


#endif
