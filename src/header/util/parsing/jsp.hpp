#ifndef __JSPREADER_HH
#define __JSPREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "format.hpp"

namespace jsp {

// template <class DT>
// DT getUb(const std::vector<DT> &duration,
//          const std::vector<std::vector<int>> &resource) {
//   DT ub{0};
//   std::vector<DT> current;
//   auto nm{resource.size()};
//   current.resize(nm + duration.size() / nm, 0);
//   std::vector<std::vector<int>> resource_of;
//   resource_of.resize(duration.size());
//   for (unsigned i = 0; i < nm; ++i) {
//     for (auto a : resource[i]) {
//       resource_of[a].push_back(i);
//     }
//   }
//   for (unsigned a{0}; a < duration.size(); ++a) {
//     resource_of[a].push_back(a / nm);
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
  auto nm{data.resources.size()};
  current.resize(nm + data.durations.size() / nm, 0);
  std::vector<std::vector<int>> resource_of;
  resource_of.resize(data.durations.size());
  for (unsigned i = 0; i < nm; ++i) {
    for (auto a : data.resources[i]) {
      resource_of[a].push_back(i);
    }
  }
  for (unsigned a{0}; a < data.durations.size(); ++a) {
    resource_of[a].push_back(a / nm);
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

ProblemInstance read_instance(const std::string &fn) {
  using std::cerr;
  ProblemInstance ret;
  ret.optimalSolution = -1;
  ret.lowerBound = 0;

  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0}, nm{0};

    int job{0};
    // int task{0};
    int dur;
    int mach;

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gotheader) {
        iss >> nj;
        iss >> nm;

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        ret.durations.reserve(nj * nm);
        ret.resources.resize(nm);
//          ret.jobs.resize(nj);

        gotheader = true;
      } else if (job < nj) {
        // job

        for (auto m = 0; m < nm; ++m) {
          iss >> mach;
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

          ret.resources[mach].push_back(ret.durations.size());
          ret.durations.push_back(dur);
        }

        ++job;
      }
    }

    for (auto job{0}; job < nj; ++job) {
        ret.constraints.push_back(
            std::make_tuple(ORIGIN, tempo::START(job * nm), 0));
      for (auto m{1}; m < nm; ++m) {
        ret.constraints.push_back(
            std::make_tuple(tempo::END(job * nm + m - 1), tempo::START(job * nm + m), 0));
      }
        ret.constraints.push_back(
            std::make_tuple(tempo::END((job+1) * nm - 1), HORIZON, 0));
    }

  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }

  return ret;
}
} // namespace jsp

#endif
