#ifndef __JSPREADER_HH
#define __JSPREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

//#include "format.hpp"

namespace jsp {



template <typename M, typename J, typename R> void parse(const std::string &fn, M &solver, J& schedule, std::vector<J>& jobs,  std::vector<R>& resources) {
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

              for (size_t m{0}; m < nm; ++m) {
                  iss >> mach;
                  iss >> dur;

                if (!iss) {
                  cerr << "ERROR: line count at line " << ln << "\n";
                  exit(1);
                }

                auto j{solver.newJob(dur, dur)};

                if (m == 0) {
                  solver.set(j.start.after(schedule.start));
                } else {
                  solver.set(jobs.back().end.before(j.start));
                  if (m == nm - 1) {
                    solver.set(j.end.before(schedule.end));
                  }
                }

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
