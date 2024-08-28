#ifndef __JSPREADER_HH
#define __JSPREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "DistanceConstraint.hpp"
#include "util/traits.hpp"

//#include "format.hpp"

namespace jsp {

template <typename M, typename J, typename R,
          tempo::concepts::same_template<tempo::DistanceConstraint> D = tempo::DistanceConstraint<int>>
void parse(const std::string &fn, M &solver, J &schedule,
           std::vector<J> &intervals, std::vector<R> &resources, std::vector<D> *precedences = nullptr) {
  using std::cerr;
  try {
      
//      std::cout << "hello\n";
      
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

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
      } else if (intervals.size() < nj * nm) {
          
//          std::cout << "job line\n";

        for (size_t m{0}; m < nm; ++m) {
          iss >> mach;
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

//          auto j{solver.newInterval(dur, dur)};
//            auto s{solver.newNumeric()};
//            auto j{solver.continuefor(solver.newNumeric(), dur)};
            
            auto s{solver.newNumeric()};
            auto j{solver.between(s, s+dur)};
            
          if (m == 0) {
            solver.set(j.start.after(schedule.start));
            if (nullptr != precedences) {
                precedences->emplace_back(j.start.id(), schedule.start.id(), 0);
            }
          } else {
            solver.set(intervals.back().end.before(j.start));
            if (nullptr != precedences) {
                precedences->emplace_back(j.start.id(), intervals.back().end.id(), 0);
            }
            if (m == nm - 1) {
              solver.set(j.end.before(schedule.end));
              if (nullptr != precedences) {
                  precedences->emplace_back(schedule.end.id(), j.end.id(), 0);
              }
            }
          }

            resources[mach].push_back(intervals.size());
          intervals.push_back(j);
//          resources[mach].push_back(j);
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
