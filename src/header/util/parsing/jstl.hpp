#ifndef __JSTLREADER_HH
#define __JSTLREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

namespace jstl {

//template <class DT> DT getUb(const ProblemInstance &data) {
//  DT ub{0};
//  for (auto d : data.durations)
//    ub += d;
//  return ub;
//}

template <typename M, typename J, typename R,
        tempo::concepts::same_template<tempo::DistanceConstraint> D = tempo::DistanceConstraint<int>>
void parse(const std::string &fn, M &solver, J &schedule,
           std::vector<J> &Intervals, std::vector<R> &resources, std::vector<D> *precedences = nullptr) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};
    bool gotsecondheader{false};
    bool stupid{true};

    std::string lt;
    std::string dump;
    int ln{1};
    size_t nj{0}, nm{0};
    int dur;
    int mach;

    int min_tl{0};
    int max_tl{0};

    int optimalSolution{-1};
    //        int lowerBound{0};

    int k{0};

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

//                        std::cout << nj << " " << nm << std::endl;

      } else if (Intervals.size() < nj * nm) {

        for (size_t m{0}; m < nm; ++m) {
          iss >> mach;
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

          auto j{solver.newInterval(dur, dur)};

          //                    std::cout << j << std::endl;

            resources[mach].push_back(Intervals.size());
          Intervals.push_back(j);
          //resources[mach].push_back(j);
        }
      } else if (not gotsecondheader and Intervals.size() == nj * nm) {
        iss >> lt;
        assert(lt == "RES=");

        iss >> optimalSolution;

        gotsecondheader = true;
      } else {

        for (size_t m{0}; m < nm; ++m) {
          if (stupid) {
            iss >> dump;
            assert(dump == "TL=0");
            min_tl = 0;
            stupid = false;
          } else {
            iss >> min_tl;
          }
          iss >> max_tl;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

          if (m == 0) {
            auto l{Intervals[k].start.after(schedule.start)};
//                                    std::cout << l << std::endl;
            solver.set(l);
            if (nullptr != precedences) {
                precedences->emplace_back(l);
            }
          } else {
            auto l{Intervals[k - 1].end.before(Intervals[k].start, min_tl)};
//                                    std::cout << l << std::endl;
            solver.set(l);
            if (nullptr != precedences) {
                precedences->emplace_back(l);
            }
            if (m == nm - 1) {
              l = Intervals[k].end.before(schedule.end);
//                                          std::cout << l << std::endl;
              solver.set(l);
              if (nullptr != precedences) {
                  precedences->emplace_back(l);
              }
            }
          }
          if (m < (nm - 1)) {
            auto l{Intervals[k + 1].start.before(Intervals[k].end, -max_tl)};
//                                    std::cout << l << std::endl;
            solver.set(l);
            if (nullptr != precedences) {
                precedences->emplace_back(l);
            }
          }
          ++k;
        }
      }
    }
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

} // namespace jstl

#endif
