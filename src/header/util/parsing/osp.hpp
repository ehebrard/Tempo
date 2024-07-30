#ifndef __OSPREADER_HH
#define __OSPREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace osp {

template <typename M, typename J, typename R>
int parse(const std::string &fn, M &model, J &schedule,
           std::vector<J> &Intervals, std::vector<R> &resources) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gothints{false};
    bool gotheader{false};

    int lowerBound;
    int optimalSolution;

    std::string lt;
    int ln{1};
    size_t nj{0}, nm{0};
    int dur;

    size_t j{0};

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gothints) {
        gothints = true;
        char _;
        iss >> lowerBound >> _ >> optimalSolution;
        continue;
      } else if (!gotheader) {
        iss >> nj;
        iss >> nm;

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        resources.resize(nj + nm);

        gotheader = true;
      } else if (j < nj) {

        for (size_t mach{0}; mach < nm; ++mach) {
          iss >> dur;

          if (!iss) {
            cerr << "ERROR: line count at line " << ln << "\n";
            exit(1);
          }

            auto s{model.newNumeric()};
            auto t{model.between(s, s+dur)};

          model.post(t.start.after(schedule.start));
          model.post(t.end.before(schedule.end));

            resources[j].push_back(Intervals.size());
            resources[nj + mach].push_back(Intervals.size());
          Intervals.push_back(t);
//          resources[j].push_back(t);
//          resources[nj + mach].push_back(t);
        }

        ++j;
      }
    }

    return optimalSolution;
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

} // namespace osp

#endif
