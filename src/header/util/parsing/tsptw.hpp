#ifndef __TSPTWREADER_HH
#define __TSPTWREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace tsptw {

template <typename M, typename J, typename R, typename D>
void parse(const std::string &fn, M &model, J &schedule,
           std::vector<J> &intervals, std::vector<int> &weights,
           std::vector<R> &resources, std::vector<D> &distances,
           const bool optional = false) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");

    bool gotheader{false};

    std::string lt; // line type
    int ln{1};
    int nj{0};

    int job{0};
    int row{0};
    int release;
    int duedate;
    int dur;

    //      std::cout << "hello\n";

    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      if (!gotheader) {
        iss >> nj;

        if (!iss) {
          cerr << "ERROR: could not parse header at line " << ln << "\n";
          exit(1);
        }

        //          std::cout << nj << " jobs\n";

        resources.resize(1);
        distances.resize(1);

        distances[0].resize(nj);
        for (auto &d : distances[0])
          d.resize(nj, 0);

        if (optional) {
          for (auto i{0}; i < nj; ++i) {
            int w;
            iss >> w;
            weights.push_back(w);
          }
        }

        gotheader = true;
      } else if (job < nj) {
        // job

        iss >> release;
        iss >> duedate;
        iss >> dur;

        if (!iss) {
          cerr << "ERROR: line count at line " << ln << "\n";
          exit(1);
        }

        auto s{model.newNumeric(release, duedate - dur)};

        tempo::Interval<int> t;
        if (optional)
          t = model.maybe_between(s, s + dur);
        else
          t = model.between(s, s + dur);

        model.post(t.start.after(schedule.start));
        model.post(t.end.before(schedule.end));

          resources[0].push_back(intervals.size());
        intervals.push_back(t);
//        resources[0].push_back(t);

        //          std::cout << " job" << job << ": " << t << std::endl;

        ++job;
      } else {

        //          std::cout << " row" << row << ":" ;

        for (auto j{0}; j < nj; ++j) {
          iss >> distances[0][row][j];

          //            std::cout << " " << distances[0][row][j];
        }

//        std::cout << std::endl;

        ++row;
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
