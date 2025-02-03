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

//        auto json{nlohmann::json::parse(file).get()};
        
        std::cout << "hello\n";
        
        
    } catch (std::exception &e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
    
} // namespace tsptw

#endif
