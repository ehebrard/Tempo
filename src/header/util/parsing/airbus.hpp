#ifndef __AIRBUSREADER_HH
#define __AIRBUSREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <nlohmann/json.hpp>

namespace airbus {


template <typename M>
void parse(const std::string &fn, M &model) {
    
    using std::cerr;
    try {
        
        std::ifstream file(fn);
        if (not file.is_open()) {
            std::stringstream ss;
            ss << "unable to open file '" << fn << "'";
            throw std::runtime_error(ss.str());
        }

        auto json{nlohmann::json::parse(file)};
        
//        std::cout << json["teams"].size() << std::endl;
        
        
//        auto origin = model.newConstant(0);
//        auto makespan = model.newNumeric(0, tempo::Constant::Infinity<int>);
        
        // from task keys to variable ids (0 is for constants, 1 is for the schedule duration/end)
//        std::map<std::string, typename tempo::SchedulingInstance<int>::ModNum> var_map;
//        auto t{2};
//        for(auto task_key : json["tasks"]) {
//            var_map[task_key] = t++;
//        }
        
//        auto schedule = model.between(origin, makespan);
        
//        auto schedule
        
        auto schedule{model.newInterval(0, tempo::Constant::Infinity<int>, 0, 0, 0,
                                    tempo::Constant::Infinity<int>)};
        
        for(auto t : json["teams"]) {
            model.post(NoOverlap(schedule));
        }
        
        
//        auto num_resources{json["teams"].size()};
//        auto num_intervals{json["tasks"].size()};
//        
//        model.declareDisjunctiveResources(num_resources);
//        model.declareCardinalityConstraints(num_intervals);
//        
//        
//        
//        auto t{0};
//        for(auto& task_key : json["tasks"]) {
//            
//            
//            auto lb{static_cast<int>(json["start_window"][task_key][0])};
//            auto ub{static_cast<int>(json["start_window"][task_key][1])};
//            auto dur{static_cast<int>(json["tasks_data"][task_key]["duration"])};
//            auto s{model.addNumeric(lb,ub)};
//            auto p{model.addConstant(dur)};
//            auto e{model.addView(s,dur)};
//            var_map[task_key] = s;
//            
//            std::cout << "task ";
//            
//            
//            var_map[task_key].display(std::cout);
//            
//            std::cout << ": dur="
//            << dur << " start window = " << lb << ".." << ub << std::endl;
//            
//            for(auto& team_key : json["teams"]) {
//                auto r{static_cast<int>(json["teams_to_index"][team_key])};
////                for(auto& task_key : json["tasks"]) {
//                for(size_t x{0}; x<num_intervals; ++x) {
//                    auto b{model.addBoolean()};
//                    auto j{model.addInterval(s,e,p,b)};
//                    
//                    model.addDisjunctiveResourceUsage(j, r);
//                    model.addCardinalityArgument(b, t);
//                }
//            }
//            
//            ++t;
//            
////            std::cout << s << " " << j << " " << var_map[task_key] << std::endl;
//        }
//        
//        for(auto& task_key : json["tasks"]) {
//            auto dur{static_cast<int>(json["tasks_data"][task_key]["duration"])};
//            std::cout << "*prec ";
//            var_map[task_key].display(std::cout);
//            std::cout << " << ";
//            model.makespan.display(std::cout);
//            std::cout << std::endl;
//            
//            model.addPrecedence(var_map[task_key], model.makespan, dur);
//            
//            for(auto& successor_key : json["successors"][task_key]) {
//                std::cout << "+prec " ;
//                
//                var_map[task_key].display(std::cout);
//                
//                std::cout << " << " ;
//                var_map[successor_key].display(std::cout);
//                
//                std::cout << std::endl;
//                
//                model.addPrecedence(var_map[task_key], var_map[successor_key], dur);
//            }
//        }
//        
//        
////        for(auto& team_key : json["teams"]) {
////            auto r{static_cast<int>(json["teams_to_index"][team_key])};
////            for(auto& task_key : json["tasks"]) {
////                model.addDisjunctiveResourceUsage(var_map[task_key], r);
////            }
////        }
//        
//        
//        for(auto &clique : json["same_allocation"]) {
//            for(auto i : clique)
//                for(auto j : clique)
//                    if(var_map[i].id < var_map[j].id) {
//                        for(size_t k{0}; k<num_resources; ++k) {
//                            model.addSameAllocation(k * num_intervals + var_map[i].id-1, k * num_intervals + var_map[j].id-1);
//                        }
//                    }
//        }
//        
//        
////        for
//        
//        
////        auto zero{model.addNumeric(0,0)};
////        model.declareDisjunctiveResources(json["teams"].size());
////        
////        for(auto task_key : json["tasks"]) {
////            auto s{model.addNumeric()};
////            auto j{model.addFixedDurationIntervalFrom(s, json["tasks_data"][task_key]["duration"])};
////        }
        
        
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
    
} // namespace tsptw

#endif
