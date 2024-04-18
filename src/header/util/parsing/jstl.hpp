#ifndef __JSTLREADER_HH
#define __JSTLREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "format.hpp"

namespace jstl {

template <class DT> DT getUb(const ProblemInstance &data) {
  DT ub{0};
  for (auto d : data.durations)
    ub += d;
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
    bool gotsecondheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0}, nm{0};

    int job{0};
    // int task{0};
    int dur;
    int mach;
    int min_tl;
    int max_tl;

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      if (line.compare(0, 3, "TL=") == 0)
        line = line.substr(3);
      std::istringstream iss(line);

      if (not gotheader) {
        iss >> nj;
        iss >> nm;

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        ret.durations.reserve(nj * nm);
        ret.resources.resize(nm);

        gotheader = true;
      } else if (not gotsecondheader and job < nj) {
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
      } else if (job == nj) {
        iss >> lt;
        assert(lt == "RES=");

        iss >> ret.optimalSolution;

        job = 0;
        gotsecondheader = true;
      } else {
        for (auto m = 0; m < nm; ++m) {
          iss >> min_tl;
          iss >> max_tl;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

          if (m < (nm - 1)) {
            ret.constraints.push_back(
                std::make_tuple(tempo::END(job * nm + m),
                                tempo::START(job * nm + m + 1), min_tl));
          } else {
            ret.constraints.push_back(
                std::make_tuple(tempo::END((job + 1) * nm - 1), HORIZON, 0));
          }

          if (m == 0) {
            ret.constraints.push_back(
                std::make_tuple(ORIGIN, tempo::START(job * nm), 0));

          } else {

            ret.constraints.push_back(
                std::make_tuple(tempo::START(job * nm + m),
                                tempo::END(job * nm + m - 1), -max_tl));
          }
        }
        ++job;
      }
    }

    //    for (auto job{0}; job < nj; ++job) {
    //        ret.constraints.push_back(
    //            std::make_tuple(ORIGIN, tempo::START(job * nm), 0));
    //      for (auto m{1}; m < nm; ++m) {
    //        ret.constraints.push_back(
    //            std::make_tuple(tempo::END(job * nm + m - 1), tempo::START(job
    //            * nm + m), 0));
    //      }
    //        ret.constraints.push_back(
    //            std::make_tuple(tempo::END((job+1) * nm - 1), HORIZON, 0));
    //    }

  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }

  return ret;
}
} // namespace jstl

#endif
