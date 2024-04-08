#ifndef __OSPREADER_HH
#define __OSPREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "format.hpp"

namespace osp {

// template <class DT>
// DT getUb(const std::vector<DT> &duration,
//          const std::vector<std::vector<int>> &resource) {
//   DT ub{0};
//   std::vector<DT> current;
//   current.resize(resource.size(), 0);
//   std::vector<std::vector<int>> resource_of;
//   resource_of.resize(duration.size());
//   for (unsigned i = 0; i < resource.size(); ++i) {
//     for (auto a : resource[i]) {
//       resource_of[a].push_back(i);
//     }
//   }
//   for (unsigned i = 0; i < duration.size(); ++i) {
//     auto d = duration[i];
//     auto first = 0;
//
//     for (auto j : resource_of[i]) {
//       if (current[j] > first)
//         first = current[j];
//     }
//
//     for (auto j : resource_of[i]) {
//       current[j] = first + d;
//       if (ub < current[j])
//         ub = current[j];
//     }
//   }
//
//   return ub;
// }

template <class DT> DT getUb(const ProblemInstance &data) {
  DT ub{0};
  std::vector<DT> current;
  current.resize(data.resources.size(), 0);
  std::vector<std::vector<int>> resource_of;
  resource_of.resize(data.durations.size());
  for (unsigned i = 0; i < data.resources.size(); ++i) {
    for (auto a : data.resources[i]) {
      resource_of[a].push_back(i);
    }
  }
  for (unsigned i = 0; i < data.durations.size(); ++i) {
    auto d = data.durations[i];
    auto first = 0;

    for (auto j : resource_of[i]) {
      if (current[j] > first)
        first = current[j];
    }

    for (auto j : resource_of[i]) {
      current[j] = first + d;
      if (ub < current[j])
        ub = current[j];
    }
  }

  return ub;
}

const ProblemInstance read_instance(const std::string &fn) {
  using std::cerr;
  ProblemInstance ret;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gothints{false};
    bool gotheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0}, nm{0};

    int job{0};
    int dur;

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gothints) {
        gothints = true;
        char _;
        iss >> ret.lowerBound  >> _ >> ret.optimalSolution;
        continue;
      } else if (!gotheader) {
        iss >> nj;
        iss >> nm;

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        ret.durations.reserve(nj * nm);
        ret.resources.resize(nj + nm);

        gotheader = true;
      } else if (job < nj) {
        // job

        for (auto mach = 0; mach < nm; ++mach) {
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

          ret.resources[job].push_back(ret.durations.size());
          ret.resources[nj + mach].push_back(ret.durations.size());
          ret.durations.push_back(dur);
        }

        ++job;
      }
    }
      
      for(unsigned i{0}; i<ret.durations.size(); ++i) {
          ret.constraints.push_back(
                                    std::make_tuple(ORIGIN, tempo::START(i), 0));
          
          ret.constraints.push_back(
                                    std::make_tuple(tempo::END(i), HORIZON, 0));
      }
      
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }

  return ret;
}

} // namespace osp

#endif
