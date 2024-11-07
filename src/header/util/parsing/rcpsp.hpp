#ifndef __RCPSPREADER_HH
#define __RCPSPREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace rcpsp {

using namespace tempo;

inline void parse(const std::string &fn, Solver<int> &model, Interval<int> &schedule,
           std::vector<Interval<int>> &intervals,
           std::vector<std::vector<size_t>> &resources,
           std::vector<std::vector<int>> &demands, std::vector<int> &capacities,
           std::vector<std::pair<int, int>> &precedences, std::vector<DistanceConstraint<int>> *distConstraints = nullptr) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    std::string lt;
    int nj{0}, nm{0};
    int dur;

    int succ;

    int dem;
    int nsucc;
    int capa;

    ifs >> nj;
    ifs >> nm;

    resources.resize(nj - 2);
    demands.resize(nj - 2);

    for (auto m{0}; m < nm; ++m) {
      ifs >> capa;
      capacities.push_back(capa);
    }

    // dummy starting job
    for (auto j{0}; j < nj; ++j) {
      ifs >> dur;
      for (auto m{0}; m < nm; ++m) {
        ifs >> dem;

        assert(j > 0 or dem == 0);

        if (j > 0 and dem > 0) {
          resources[j - 1].push_back(m);
          demands[j - 1].push_back(dem);
        }
      }
      ifs >> nsucc;
      for (auto i{0}; i < nsucc; ++i) {
        ifs >> succ;
        precedences.emplace_back(j - 1, succ - 2);
      }

      if (j > 0 and j < nj - 1) {
        auto s{model.newNumeric()};
        auto t{model.between(s, s + dur)};
        intervals.push_back(t);
      }
    }

    if (nullptr != distConstraints) {
        distConstraints->reserve(precedences.size());
    }

    for (auto prec : precedences) {
      auto x{prec.first};
      auto y{prec.second};
      DistanceConstraint<int> dc;

//      std::cout << x << " < " << y << " (" << nj - 2 << ")\n";

      if (x >= nj - 2) {
        continue;
      } else if (x < 0) {
        dc = schedule.start.before(intervals[y].start);
        model.post(dc);
        if (nullptr != distConstraints) {
            distConstraints->emplace_back(dc);
        }
      } else if (y >= nj - 2) {
        dc = intervals[x].end.before(schedule.end);
        model.post(dc);
          if (nullptr != distConstraints) {
              distConstraints->emplace_back(dc);
          }
      } else {
        dc = intervals[x].end.before(intervals[y].start);
        model.post(dc);
          if (nullptr != distConstraints) {
              distConstraints->emplace_back(dc);
          }
      }

//      std::cout << dc << std::endl;
    }

//    exit(1);

  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

} // namespace osp

#endif
