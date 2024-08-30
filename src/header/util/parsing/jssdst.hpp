#ifndef __JSSDSTREADER_HH
#define __JSSDSTREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace jssdst {

template <typename M, typename J, typename R, typename D>
void parse(const std::string &fn, M &model, J &schedule,
           std::vector<J> &intervals,
           std::vector<R> &resources, std::vector<D> &distances) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0};
      int nm{0};
      int nf{0};
      int n{0};

    int job{0};
      int row{0};
    int dur;
      int mach;
      int fam;
      
      std::vector<int> family;
      std::vector<std::vector<int>> distance_matrix;

//          std::cout << "hello\n";

      for (std::string line; getline(ifs, line); ++ln) {
          if (line[0] == '#')
              continue;
          
          std::istringstream iss(line);
          
          if (!gotheader) {
              iss >> nm;
              iss >> nj;
              iss >> nf;
              
              if (!iss) {
                  cerr << "ERROR: could not parse header at line " << ln << "\n";
                  exit(1);
              }
              
              //                  std::cout << nj << " jobs " << nm << " machines " << nf << " families \n";
              
              resources.resize(nm);
              distances.resize(nm);
              
              //          family.resize(nm*nj);
              distance_matrix.resize(nf);
              for (auto & row : distance_matrix)
                  row.resize(nf, 0);
              
              
              for(auto& matrix : distances) {
                  matrix.resize(nj);
                  for (auto & row : matrix)
                      row.resize(nj, 0);
              }
              
              gotheader = true;
          } else if (job < nj) {
              // job
              
              iss >> n;
              
              //          std::cout << "job " << job << ": length=" << n << std::endl;
              
              for(int i{0}; i<n; ++i) {
                  
                  iss >> dur;
                  iss >> mach;
                  iss >> fam;
                  
                  //              std::cout << " (" << dur << "|" << mach << "|" << fam << ")";
                  
                  if (!iss) {
                      cerr << "ERROR: line count at line " << ln << "\n";
                      exit(1);
                  }
                  
                  int task{static_cast<int>(intervals.size())};
                  resources[mach-1].push_back(task);
                  family.push_back(fam-1);
                  
                  auto s{model.newNumeric()};
                  auto j{model.between(s, s + dur)};
                  
                  if (i == 0) {
                      model.set(j.start.after(schedule.start));
                      //                if (nullptr != precedences) {
                      //                    precedences->emplace_back(j.start.id(), schedule.start.id(), 0);
                      //                }
                  } else {
                      model.set(intervals.back().end.before(j.start));
                      //                if (nullptr != precedences) {
                      //                    precedences->emplace_back(j.start.id(), intervals.back().end.id(), 0);
                      //                }
                      if (i == nm - 1) {
                          model.set(j.end.before(schedule.end));
                          //                  if (nullptr != precedences) {
                          //                      precedences->emplace_back(schedule.end.id(), j.end.id(), 0);
                          //                  }
                      }
                  }
                  
                  intervals.push_back(j);
              }
              ++job;
          } else if(row < nf) {
              
              //          for (auto i{0}; i < nf; ++i) {
              for (auto j{0}; j < nf; ++j) {
                  iss >> distance_matrix[row][j];
//                  std::cout << " " << distance_matrix[row][j] ;
              }
//              std::cout << std::endl;
              ++row;
              //          }
              //          std::cout << std::endl;
          }
      }
          
          int k{0};
          for(auto & tasks : resources) {
              
//              for(auto t : tasks)
//                  std::cout << " " << t;
//              std::cout << std::endl;
              
              assert(static_cast<int>(tasks.size()) == nj);
              
//              for(size_t i{0}; i<tasks.size(); ++i) {
//                  std::cout << tasks[i] << ": " << family[tasks[i]] << std::endl;
//              }
              
              for(size_t i{0}; i<tasks.size(); ++i) {
                  for(size_t j{0}; j<tasks.size(); ++j) {
                      if(i != j) {
                          distances[k][i][j] = distance_matrix[family[tasks[i]]][family[tasks[j]]];
                      }
//                      std::cout << " " << distances[k][i][j] ;
                  }
//                  std::cout << std::endl;
              }
//              std::cout << std::endl;
              
              ++k;
          }

//          exit(1);
//      }
//    }

//        std::cout << model << std::endl;
//      exit(1);

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
