#ifndef __OSPREADER_HH
#define __OSPREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace osp {

template <typename M, typename J, typename R>
void parse(const std::string &fn, M &model, J &schedule,
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

    //        std::cout << "makespan\n";
    //        auto schedule{model.newInterval()};
    //        schedule = model.newInterval();
    model.set(schedule.start.after(0));
    model.set(schedule.start.before(0));

    size_t j{0};

    //        int Interval{0};
    for (std::string line; getline(ifs, line); ++ln) {
      if (line[0] == '#')
        continue;

      std::istringstream iss(line);

      //            std::cout << "line " << ln << " " << Intervals.size() <<
      //            "Intervals\n";

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

          //                std::cout << "new task (" << dur << ")\n";
          auto t{model.newInterval(dur, dur)};

          model.set(t.start.after(schedule.start));
          model.set(t.end.before(schedule.end));

          Intervals.push_back(t);
          resources[j].push_back(t);
          resources[nj + mach].push_back(t);
        }

        ++j;
      }
    }
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

//
// template <typename S, typename M>
// void parse(const std::string &fn, S &solver, M &model) {
//  using std::cerr;
//  try {
//    std::ifstream ifs(fn);
//    if (!ifs)
//      throw std::runtime_error("Could not open file for reading");
//
//    bool gothints{false};
//    bool gotheader{false};
//
//    int lowerBound;
//    int optimalSolution;
//
//    std::string lt;
//    int ln{1};
//    size_t nj{0}, nm{0};
//    int dur;
//
//    //        std::cout << "makespan\n";
//    //        auto schedule{model.newInterval()};
//    //        schedule = model.newInterval();
//
//      auto l{schedule.start.after(0)};
//      model.precedences.push_back(l);
//      solver.set(l);
//
//      auto l{schedule.start.before(0)};
//      model.precedences.push_back(l);
//      solver.set(l);
//
//
//    size_t j{0};
//
//    //        int Interval{0};
//    for (std::string line; getline(ifs, line); ++ln) {
//      if (line[0] == '#')
//        continue;
//
//      std::istringstream iss(line);
//
//      //            std::cout << "line " << ln << " " << Intervals.size() <<
//      //            "Intervals\n";
//
//      if (!gothints) {
//        gothints = true;
//        char _;
//        iss >> lowerBound >> _ >> optimalSolution;
//        continue;
//      } else if (!gotheader) {
//        iss >> nj;
//        iss >> nm;
//
//        if (!iss) {
//          cerr << "ERROR: could not parse header at line " << ln << "\n";
//          exit(1);
//        }
//
//        resources.resize(nj + nm);
//
//        gotheader = true;
//      } else if (j < nj) {
//
//        for (size_t mach{0}; mach < nm; ++mach) {
//          iss >> dur;
//
//          if (!iss) {
//            cerr << "ERROR: line count at line " << ln << "\n";
//            exit(1);
//          }
//
//          //                std::cout << "new task (" << dur << ")\n";
////          auto t{model.newInterval(dur, dur)};
////
////          model.set(t.start.after(schedule.start));
////          model.set(t.end.before(schedule.end));
//
//            auto t{solver.newInterval(dur, dur)};
//            model.push_back(t);
//
//            auto l{t.start.after(schedule.start)};
//            model.push_back(l)
//
//            solver.set(l);
//
//            l = t.end.before(schedule.end);
//            solver.set(t);
//
//
//          Intervals.push_back(t);
//          resources[j].push_back(t);
//          resources[nj + mach].push_back(t);
//        }
//
//        ++j;
//      }
//    }
//  } catch (std::exception &e) {
//    std::cout.flush();
//    cerr << "ERROR: " << e.what() << std::endl;
//    exit(1);
//  }
//}

} // namespace osp

#endif
