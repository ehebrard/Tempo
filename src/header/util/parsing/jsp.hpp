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
            auto d = j.start.after(schedule.start);
            solver.set(d);
            if (nullptr != precedences) {
                precedences->emplace_back(d);
            }
          } else {
            auto d = intervals.back().end.before(j.start);
            solver.set(d);
            if (nullptr != precedences) {
                precedences->emplace_back(d);
            }
            if (m == nm - 1) {
              auto d = j.end.before(schedule.end);
              solver.set(d);
              if (nullptr != precedences) {
                  precedences->emplace_back(d);
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



template <typename M>
void parse(const std::string &fn, M &model) {
  using std::cerr;
  try {
      
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt;
    int ln{1};
    size_t nj{0}, nm{0};
    int dur;
    int mach;
      
      size_t num_intervals{0};
      
      
      auto zero{model.addNumeric(0,0)};
//      auto horizon{model.addNumeric(0)};
//      model.addInterval(zero, horizon);

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

        model.declareDisjunctiveResources(nm);

        gotheader = true;
      } else if (num_intervals < nj * nm) {
  
          auto prev{zero};
          auto prev_dur{0};
        for (size_t m{0}; m < nm; ++m) {
          iss >> mach;
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }
  
            auto s{model.addNumeric()};
            auto j{model.addFixedDurationIntervalFrom(s, dur)};
            
            model.addPrecedence(prev, s, prev_dur);
            
            prev = s;
            prev_dur = dur;
  
            model.addDisjunctiveResourceUsage(j, mach);
        }
//          model.addPrecedence(prev, horizon, prev_dur);
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
