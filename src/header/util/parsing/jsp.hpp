#ifndef __JSPREADER_HH
#define __JSPREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

//#include "format.hpp"

namespace jsp {



//template <class DT> DT getUb(const ProblemInstance &data) {
//  DT ub{0};
//  std::vector<DT> current;
//  auto nm{data.resources.size()};
//  auto nj{data.durations.size() / nm};
//  current.resize(nm + nj, 0);
//  std::vector<std::vector<int>> resource_of(data.durations.size());
//  for (unsigned i = 0; i < nm; ++i) {
//    for (auto a : data.resources[i]) {
//      resource_of[a].push_back(i);
//    }
//  }
//  for (unsigned a{0}; a < data.durations.size(); ++a) {
//    resource_of[a].push_back(a / nm);
//  }
//
//  std::vector<DT> trail(data.durations.size());
//  for (size_t j{0}; j < nj; ++j) {
//    trail[j * nm + nm - 1] = 0;
//  }
//  for (auto i{nm}; --i > 0;) {
//    for (size_t j{0}; j < nj; ++j) {
//      trail[j * nm + i - 1] = trail[j * nm + i] + data.durations[j * nm + i];
//    }
//  }
//
//  std::vector<int> insertion_order(nm * nj);
//  for (size_t i{0}; i < nm * nj; ++i) {
//    insertion_order[i] = i;
//  }
//
//  std::sort(insertion_order.begin(), insertion_order.end(),
//            [&](const int x, const int y) {
//              return (((x % nm) < (y % nm)) or
//                      ((x % nm) == (y % nm) and trail[x] > trail[y]));
//            });
//
//  //        // proceed by job rank
//  //    for(unsigned i{0}; i<nm; ++i) {
//  //
//  //    }
//  //
//
//  for (unsigned z = 0; z < data.durations.size(); ++z) {
//    auto i{insertion_order[z]};
//    auto d = data.durations[i];
//    auto first = 0;
//
//    for (auto j : resource_of[i]) {
//      if (current[j] > first)
//        first = current[j];
//    }
//
//    for (auto j : resource_of[i]) {
//      current[j] = first + d;
//      if (ub < current[j])
//        ub = current[j];
//    }
//  }
//
//  return ub;
//}

//ProblemInstance read_instance(const std::string &fn) {
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
//    bool gotheader{false};
//
//    std::string lt; // line type
//    int ln{1};
//    int nj{0}, nm{0};
//
//    int job{0};
//    // int task{0};
//    int dur;
//    int mach;
//
//    for (std::string line; getline(ifs, line); ++ln) {
//      if (line[0] == '#')
//        continue;
//
//      std::istringstream iss(line);
//
//      if (!gotheader) {
//        iss >> nj;
//        iss >> nm;
//
//        if (!iss) {
//          cerr << "ERROR: could not parse header at line " << ln << "\n";
//          exit(1);
//        }
//
//        ret.durations.reserve(nj * nm);
//        ret.resources.resize(nm);
////          ret.jobs.resize(nj);
//
//        gotheader = true;
//      } else if (job < nj) {
//        // job
//
//        for (auto m = 0; m < nm; ++m) {
//          iss >> mach;
//          iss >> dur;
//
//          if (!iss) {
//            cerr << "ERROR: line count at line " << ln << "\n";
//            exit(1);
//          }
//
//          ret.resources[mach].push_back(ret.durations.size());
//          ret.durations.push_back(dur);
//        }
//
//        ++job;
//      }
//    }
//
//    for (auto job{0}; job < nj; ++job) {
//        ret.constraints.push_back(
//            std::make_tuple(ORIGIN, tempo::START(job * nm), 0));
//      for (auto m{1}; m < nm; ++m) {
//        ret.constraints.push_back(
//            std::make_tuple(tempo::END(job * nm + m - 1), tempo::START(job * nm + m), 0));
//      }
//        ret.constraints.push_back(
//            std::make_tuple(tempo::END((job+1) * nm - 1), HORIZON, 0));
//    }
//
//  } catch (std::exception &e) {
//    std::cout.flush();
//    cerr << "ERROR: " << e.what() << std::endl;
//    exit(1);
//  }
//
//  return ret;
//}



template <typename M, typename J, typename R> void parse(const std::string &fn, M &solver, J& schedule, std::vector<J>& jobs,  std::vector<R>& resources) {
    using std::cerr;
    try {
        std::ifstream ifs(fn);
        if (!ifs)
            throw std::runtime_error("Could not open file for reading");
        
//        bool gothints{false};
        bool gotheader{false};
        
//        int lowerBound;
//        int optimalSolution;
        
        std::string lt;
        int ln{1};
        size_t nj{0}, nm{0};
        int dur;
        int mach;

        
        solver.set(schedule.start.after(0));
        solver.set(schedule.start.before(0));
        

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
                                
                resources.resize(nm);
                
                gotheader = true;
            } else if (jobs.size() < nj * nm) {
              //                J prev;

              for (size_t m{0}; m < nm; ++m) {
                  iss >> mach;
                  iss >> dur;

                if (!iss) {
                  cerr << "ERROR: line count at line " << ln << "\n";
                  exit(1);
                }

//                std::cout << "new task (" << dur << ")\n";
                auto j{solver.newJob(dur, dur)};

                if (m == 0) {
                  solver.set(j.start.after(schedule.start));
                } else {
                  solver.set(jobs.back().end.before(j.start));
                  if (m == nm - 1) {
                    solver.set(j.end.before(schedule.end));
                  }
                }

                //                    prev = j;

                jobs.push_back(j);
                resources[mach].push_back(j);
              }
            }
        }
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}


} // namespace jsp

#endif
