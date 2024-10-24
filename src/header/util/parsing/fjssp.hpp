#ifndef __FJSSPREADER_HH
#define __FJSSPREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace fjssp {

void parse(const std::string &fn,
           int& num_machine,
           std::vector<int> &job_size,
           std::vector<std::vector<int>> &resource, // for each interval, the resource(s) it requires
           std::vector<std::vector<int>> &duration // for each interval, the duration on the correponding resource
           ) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0};
      int nt{0};
      int no{0};
      int nm{0};
      int mach{0};

    int job{0};
//    int row{0};
    int dur;
      int dump;
      int ti{0};

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gotheader) {
        iss >> nj;
          iss >> nm;
          iss >> dump;
          
          
          num_machine = nm;
          

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        //          std::cout << nj << " jobs\n";

//        resources.resize();
          job_size.resize(nj);

        gotheader = true;
      } else if (job < nj) {
        // job

          iss >> nt;
          
          resource.resize(resource.size() + nt);
          duration.resize(duration.size() + nt);
          job_size[job] = nt;
          
          for(auto t{0}; t<nt; ++t) {
              iss >> no;
              for(auto o{0}; o<no; ++o) {
                  iss >> mach;
                  iss >> dur;
                  resource[ti].push_back(mach-1);
                  duration[ti].push_back(dur);
//                  jobs[job].push_back(ti);
              }
              ++ti;
          }

        if (!iss) {
          cerr << "ERROR: line count at line " << ln << "\n";
          exit(1);
        }

           ++job;
      }
    }

    //      std::cout << model << std::endl;

  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

//
// ProblemInstance read_instance(const std::string &fn) {
//  using std::cerr;
//  ProblemInstance ret;
//  ret.optimalSolution = -1;
//  ret.lowerBound = 0;
//
//  try {
//    std::ifstream ifs(fn);
//    if (!ifs)
//      throw std::runtime_error("Could not open file for reading");
//
//      bool gotheader{false};
//
//      std::string lt; // line type
//      int ln{1};
//      int nj{0};
//
//      int job{0};
//      int row{0};
//      int release;
//      int duedate;
//      int dur;
//
//    for (std::string line; getline(ifs, line); ++ln) {
//      if (line[0] == '#')
//        continue;
//
//      std::istringstream iss(line);
//
//      if (!gotheader) {
//        iss >> nj;
//
//        if (!iss) {
//          cerr << "ERROR: could not parse header at line " << ln << "\n";
//          exit(1);
//        }
//
//        ret.durations.reserve(nj);
//        ret.resources.resize(1);
//        ret.resources[0].transition.resize(nj);
//
//        gotheader = true;
//      } else if (job < nj) {
//        // job
//
//        iss >> release;
//        iss >> duedate;
//        iss >> dur;
//
//        if (!iss) {
//          cerr << "ERROR: line count at line " << ln << "\n";
//          exit(1);
//        }
//
//        ret.resources[0].push_back(job);
//        ret.durations.push_back(dur);
//        ret.constraints.push_back(
//            std::make_tuple(ORIGIN, tempo::START(job), release));
//        ret.constraints.push_back(
//            std::make_tuple(tempo::END(job), ORIGIN, -duedate));
//
//        ++job;
//      } else {
//
//        ret.resources[0].transition[row].resize(nj);
//        for (auto j{0}; j < nj; ++j) {
//          iss >> ret.resources[0].transition[row][j];
//        }
//        ++row;
//      }
//    }
//  } catch (std::exception &e) {
//    std::cout.flush();
//    cerr << "ERROR: " << e.what() << std::endl;
//    exit(1);
//  }
//
//  return ret;
//}
} // namespace tsptw

#endif
