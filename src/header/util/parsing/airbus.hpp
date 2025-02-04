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
        
        
        std::map<std::string, int> task_map;
        auto t{0};
        for(auto task_key : json["tasks"]) {
            task_map[task_key] = t++;
        }
        
        
        for(auto& task_key : json["tasks"]) {
            std::cout << "task " << task_map[task_key] << ": dur=" << json["tasks_data"][task_key]["duration"] << " start window = " << static_cast<int>(json["start_window"][task_key][0]) << ".."
            << static_cast<int>(json["start_window"][task_key][1]) << std::endl;
            
            auto lb{static_cast<int>(json["start_window"][task_key][0])};
            auto ub{static_cast<int>(json["start_window"][task_key][1])};
            auto dur{static_cast<int>(json["tasks_data"][task_key]["duration"])};
            auto s{model.addNumeric(lb,ub)};
            auto j{model.addFixedDurationOptionalIntervalFrom(s,dur)};
            
            std::cout << s << " " << j << " " << task_map[task_key] << std::endl;
        }
        
        for(auto& task_key : json["tasks"]) {
            auto dur{static_cast<int>(json["tasks_data"][task_key]["duration"])};
            for(auto& successor_key : json["successors"][task_key]) {
                std::cout << "prec " << task_map[task_key] << " << " << task_map[successor_key] << std::endl;
                
                model.addPrecedence(model.getStart(task_map[task_key]), model.getStart(task_map[successor_key]), dur);
            }
        }
        
        model.declareDisjunctiveResources(json["teams"].size());
        for(auto& team_key : json["teams"]) {
            auto r{static_cast<int>(json["teams_to_index"][team_key])};
            for(auto& task_key : json["tasks"]) {
                model.addDisjunctiveResourceUsage(task_map[task_key], r);
            }
        }
        
        
//        for
        
        
//        auto zero{model.addNumeric(0,0)};
//        model.declareDisjunctiveResources(json["teams"].size());
//        
//        for(auto task_key : json["tasks"]) {
//            auto s{model.addNumeric()};
//            auto j{model.addFixedDurationIntervalFrom(s, json["tasks_data"][task_key]["duration"])};
//        }
        
        
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
    
} // namespace tsptw

#endif
