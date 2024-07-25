#ifndef __PATHREADER_HH
#define __PATHREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


#include "Constant.hpp"


namespace path {

template <typename M, typename J, typename R>
void parse(const std::string &fn, M &solver, J &schedule,
           std::vector<J> &intervals, std::vector<R> &resources) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt;
    int ln{1};
    int nj{0}, nm{0};
    int dur;
    int mach;
      int path = 0;
      int path_length;
      int nw{0};
      int np{0};
      int w;
      
      std::vector<bool> nowait;


    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gotheader) {
        iss >> nj;
        iss >> nm;
          
//          std::cout << nj << " jobs / " << nm << " machines\n";
          
        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }
          
          nowait.resize(nm+1, false);
          
          
          
//          for(auto i{0}; i<nm; ++i) {
//              iss >> nw;
//              nowait[i] = nw;
//              
//              if (!iss) {
//                cerr << "ERROR: could not parse header at line " << ln << "\n";
//                exit(1);
//              }
//              
////              std::cout << nowait[i] ;
//          }
////          std::cout << std::endl;
          
          iss >> nw;

          for(auto i{0}; i<nw; ++i) {
              iss >> w;
              nowait[w] = true;
              
              if (!iss) {
                cerr << "ERROR: could not parse header at line " << ln << "\n";
                exit(1);
              }
              
//              std::cout << nowait[i] ;
          }
//          std::cout << std::endl;
          

        resources.resize(nm);
//          resource_ptr.resize(nm);

        gotheader = true;
      } else if (path++ < nj) {
          
//          ++path;

          iss >> path_length;
          
//          std::cout << "path of length " << path_length << std::endl;
          
        for (auto m{0}; m < path_length; ++m) {
          iss >> mach;
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

    
//          auto j{solver.newInterval(
//              dur, (nowait[mach] ? dur : tempo::Constant::Infinity<int>))};
//          if (m == 0) {
//            solver.post(j.start.after(schedule.start));
//          } else {
//            solver.post(intervals.back().end.before(j.start));
//            solver.post(j.start.before(intervals.back().end));
//            if (m == path_length - 1) {
//              solver.set(j.end.before(schedule.end));
//            }
//          }
            
            tempo::Interval<int> j;
            if (m == 0) {
                auto s{solver.newNumeric()};
                if(nowait[mach])
                    j = solver.between(s, s+dur);
                else {
                    auto d{solver.newNumeric(dur, tempo::Constant::Infinity<int>)};
                    j = solver.continuefor(s, d);
                }
              solver.post(j.start.after(schedule.start));
            } else {
                if(nowait[mach])
                    j = solver.between(intervals.back().end, intervals.back().end+dur);
                else {
                    auto d{solver.newNumeric(dur, tempo::Constant::Infinity<int>)};
                    j = solver.continuefor(intervals.back().end, d);
                }
              if (m == path_length - 1) {
                solver.set(j.end.before(schedule.end));
              }
            }
            resources[mach].push_back(intervals.size());
          intervals.push_back(j);
//          resources[mach].push_back(j);
            
        }
      } else {
          iss >> np;
          for(auto i{0}; i<np; ++i) {
              int x, y;
              iss >> x;
              iss >> y;
//              std::cout << x << " " << y+1  << " and  " << x+1 << " " << y << std::endl;
              
              solver.addToSearch(solver.newDisjunct(intervals[x-1].start.after(intervals[y].end), intervals[y-1].start.after(intervals[x].end)));
          }
//          int count{0};
//          std::cout << "swaps!\n";
//          while (iss) {
//              int x, y;
//              iss >> x;
//              
//              if(not iss)
//                  break;
//              iss >> y;
//              if(not iss)
//                  std::cout << "not right\n";
//              else
//                  std::cout << ++count << ": " << x << " " << y+1  << " and  " << x+1 << " " << y << std::endl;
//          }
//          exit(1);
      }
    }
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

} // namespace path

#endif
