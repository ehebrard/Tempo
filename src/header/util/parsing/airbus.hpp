#ifndef __AIRBUSREADER_HH
#define __AIRBUSREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <nlohmann/json.hpp>

#include "Model.hpp"

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
        auto makespan{model.newNumeric(0, tempo::Constant::Infinity<int>)};
        auto origin{model.newConstant(0)};
        auto schedule{model.between(origin, makespan)};
        
    
        std::vector<tempo::NoOverlapExpression<int>> R;
        for(auto t : json["teams"]) {
            auto r{NoOverlap(schedule)};
            R.push_back(r);
        }
        
        std::map<std::string, size_t> interval_map;
        std::vector<tempo::Interval<int>> intervals; // the actual intervals, affectations are copies
        std::map<std::string, tempo::BooleanVar<int>> affectation_map;
        
        std::vector<std::vector<tempo::BooleanVar<int>>> alternatives(json["tasks"].size());
        
//        std::cout << "init: " << model.boolean.size() << " Boolean vars and " << model.numeric.size() << " numeric vars\n";

        for(auto& task_key : json["tasks"]) {
            
            auto
            lb_s{static_cast<int>(json["start_window"][task_key][0])};
            auto
            ub_s{static_cast<int>(json["start_window"][task_key][1])};
            auto
            lb_e{static_cast<int>(json["end_window"][task_key][0])};
            auto
            ub_e{static_cast<int>(json["end_window"][task_key][1])};
            auto
            dur{static_cast<int>(json["tasks_data"][task_key]["duration"])};
            auto s{model.newNumeric(std::max(lb_s,lb_e-dur),std::min(ub_s, ub_e-dur))};
            auto j{model.between(s, s+dur)};
            interval_map[task_key] = intervals.size();
                            
            for(auto& team_key : json["teams"]) {
                auto r{static_cast<int>(json["teams_to_index"][team_key])};
                auto oj{model.maybe_between(s, s+dur)};
                oj.require(R[r]);
                alternatives[intervals.size()].push_back(oj.exist);
            }
                                                                    
            model.post(Cardinality(alternatives[intervals.size()],1,1));
            intervals.push_back(j);
            
            
//            std::cout << "task " << task_key << ": " << model.boolean.size() << " Boolean vars and " << model.numeric.size() << " numeric vars\n";
            
        }
 
//        auto k{0};
        for(auto &r : R) {
            model.post(r);
            
//            std::cout << "resource " << k++ << ": " << model.boolean.size() << " Boolean vars and " << model.numeric.size() << " numeric vars\n";
        }
        
        
        
        for(auto& task_key : json["tasks"]) {
            model.post(intervals[interval_map[task_key]].end <= makespan);
            for(auto& successor_key : json["successors"][task_key]) {
                model.post(intervals[interval_map[task_key]].end <= intervals[interval_map[successor_key]].start);
            }
        }
        
        for(auto& tasks : json["same_allocation"]) {
            for(auto ti{tasks.begin()}; ti!=tasks.end(); ++ti) {
                auto i{interval_map[*ti]};
                for(auto tj{ti+1}; tj!=tasks.end(); ++tj) {
                    auto j{interval_map[*tj]};
                    for(size_t k{0}; k<json["teams"].size(); ++k) {
                        model.post(alternatives[i][k] == alternatives[j][k]);
                    }
                }
            }
        }
        
        for(auto& team_key : json["teams"]) {
            auto r{static_cast<int>(json["teams_to_index"][team_key])};
            auto& windows{json["calendar"][team_key]};
            for(size_t i{0}; i<intervals.size(); ++i) {
                auto affectation_ir{alternatives[i][r]};
                model.post(affectation_ir.implies(intervals[i].end >= static_cast<int>(windows[0][0])));
                for(size_t k{1}; k<windows.size(); ++k) {
                    int end_prev{windows[k-1][1]};
                    int beg_next{windows[k][0]};
                    model.post(affectation_ir.implies((intervals[i].end <= end_prev) || (intervals[i].start >= beg_next)));
                }
                model.post(affectation_ir.implies((intervals[i].end <= static_cast<int>(windows.back()[1]))));
            }
        }
        
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
    
} // namespace tsptw

#endif
