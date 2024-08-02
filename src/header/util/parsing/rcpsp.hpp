#ifndef __RCPSPREADER_HH
#define __RCPSPREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace rcpsp {

using namespace tempo;

void parse(const std::string &fn,
          Solver<int> &model,
          Interval<int> &schedule,
           std::vector<Interval<int>> &intervals,
          std::vector<std::vector<size_t>> &resources,
          std::vector<std::vector<int>> &demands,
           std::vector<int> &capacities
           ) {
    using std::cerr;
    try {
        std::ifstream ifs(fn);
        if (!ifs)
            throw std::runtime_error("Could not open file for reading");
        
        bool gotnres{false};
        bool gotnjobheader{false};
        bool gotnjob{false};
        bool gotjobcrap{false};
        bool gotrescrap{false};
        bool jobdur{false};
        int gotcapacrap{3};
    
        std::string lt;
        std::string probe;
//        int ln{1};
        int nj{0}, nm{0};
        int dur;
        
        int job{0};
        int succ;
        
        int dem;
        int nsucc;
        int capa;
        
        std::vector<int> edges;
        
        for (std::string line; getline(ifs, line);) {
            
            std::istringstream iss(line);
            
//            std::cout << "line " << ln << ": " << line << std::endl;
            
            if(not gotnres) {
//                std::cout << " ** header crap\n";
                iss >> probe;
                if(probe == "-") {
                    iss >> probe;
                    if(probe == "renewable") {
                        iss >> probe;
                        assert(probe == ":");
                        iss >> nm;
                        gotnres = true;
                    }
                }
            } else if(not gotnjobheader) {
//                std::cout << " ** njob header crap\n";
                iss >> probe;
                if(probe == "pronr.") {
                    gotnjobheader = true;
                }
            } else if(not gotnjob) {
//                std::cout << " ** njob crap\n";
                iss >> probe;
                iss >> nj;
                gotnjob = true;
                resources.resize(nj);
                demands.resize(nj);
            } else if(not gotjobcrap) {
//                std::cout << " ** job crap again\n";
                iss >> probe;
                if(probe == "jobnr.") {
                    gotjobcrap = true;
                }
            } else if(not jobdur and job < nj+2){
//                std::cout << " ** job " << job << "\n";
                iss >> job;
                iss >> probe;
                assert(probe == "1");
                iss >> nsucc;
                for(auto i{0}; i<nsucc; ++i) {
                    iss >> succ;
                    edges.push_back(job);
                    edges.push_back(succ);
                }
            } else if(not gotrescrap) {
                jobdur = true;
//                std::cout << " ** res crap\n";
                iss >> probe;
                if(probe.substr(0,3) == "---") {
                    gotrescrap = true;
                    job = 0;
                }
            } else if(job < nj+2) {
//                std::cout << " ** job/res " << job << "\n";
                iss >> job;
                iss >> probe;
                assert(probe == "1");
                iss >> dur;
                if(dur > 0) {
                    assert(static_cast<size_t>(job) == intervals.size()+2);
                    auto s{model.newNumeric()};
                    auto t{model.between(s, s+dur)};
                    intervals.push_back(t);
                }
                for(auto i{0}; i<nm; ++i) {
                    iss >> dem;
                    if(dem > 0) {
                        resources[job-2].push_back(i);
                        demands[job-2].push_back(dem);
                    }
                }
            } else if(gotcapacrap > 0) {
//                std::cout << " ** capa crap\n";
                --gotcapacrap;
            } else {
//                std::cout << " ** capacities\n";
                for(auto i{0}; i<nm; ++i) {
                    iss >> capa;
                    capacities.push_back(capa);
                }
                break;
            }
        }
                
        for(unsigned i{0}; i<edges.size(); i+=2) {
            auto x{edges[i]-2};
            auto y{edges[i+1]-2};
                       
            if(x<0) {
                model.post(schedule.start.before(intervals[y].start));
            } else if(y>=nj) {
                model.post(intervals[x].end.before(schedule.end));
            } else {
                model.post(intervals[x].end.before(intervals[y].start));
            }
        }
        
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}

} // namespace osp

#endif
