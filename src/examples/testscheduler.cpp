/*************************************************************************
minicsp

Copyright 2010--2011 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/


#include <iostream>


#include "Scheduler.hpp"


using namespace tempo;




int main(int argc, char *argv[]) {
    
    Options opt = tempo::parse(argc, argv);
    
    int n{3};
    
    Scheduler<int> S(opt);
    
    std::cout << "\nempty:\n";
    std::cout << S << std::endl;
    
    for(auto i{0}; i<n*n; ++i) {
        int p_min{static_cast<int>((tempo::random()%50)+1)};
        int p_max{static_cast<int>((tempo::random()%10)+p_min)};
        S.newTask(p_min, p_max);
        
//        std::cout << S << std::endl;
    }
    
    for(auto j{0}; j<n; ++j) {
        for(auto i{1}; i<n; ++i) {
            S.newPrecedence(END(3*j+i-1), START(3*j+i), tempo::random()%10);
            
//            std::cout << S << std::endl;
        }
        S.newPrecedence(ORIGIN, START(3*j), tempo::random()%10);
        S.newPrecedence(END(3*j+2), HORIZON, tempo::random()%10);
        
//        std::cout << S << std::endl;
    }
    
    S.setUpperBound(1000);
    
    std::cout << "\nprecedences:\n";
    std::cout << S << std::endl;
    
    
    for(auto j{0}; j<n; ++j) {
        for(auto i{0}; i<n; ++i) {
            for(auto k{i+1}; k<n; ++k) {
                S.newDisjunct(START(3*i+j), END(3*k+j), 0, START(3*k+j), END(3*i+j), 0);
//                S.newVariable(START(3*i+j), END(3*k+j), 0);
//                S.newVariable(START(3*k+j), END(3*i+j), 0);
            }
        }
    }
    
    std::cout << "\ndisjuncts (save):\n";
    std::cout << S << std::endl;
    
    
    S.saveState();
    
    // e0 < s3
    S.set(1,true);
        
    // not e1 < s7 -> s7 < e1
    S.set(9,false);
        
    // e3 < s6
    S.set(5,true);
    
    // e8 < s2
    S.set(14,true);
    
    // e5 < s8
    S.set(17,true);
    
    std::cout << "\ndecisions:\n";
    std::cout << S << std::endl;
    
    
    S.propagate();
    
    std::cout << "\npropag:\n";
    std::cout << S << std::endl;
    
    S.restoreState(0);
    
    std::cout << "\nrestore:\n";
    std::cout << S << std::endl;
    
    
    try {
        S.setUpperBound(200);
        
//        std::cout << S << std::endl;
        
        S.propagate();
        
        std::cout << "\ntighter (save):\n";
        std::cout << S << std::endl;
        
    } catch(std::exception& e) {
        std::cout << e.what() << std::endl;
        exit(1);
    }
    S.saveState();
    
    // e0 < s3
    S.set(1,true);
        
    // not e1 < s7 -> s7 < e1
    S.set(9,false);
        
    // e3 < s6
    S.set(5,true);
    
    // e8 < s2
    S.set(14,true);
    
    // e5 < s8
    S.set(17,true);
    
    std::cout << "\ndecisions:\n";
    std::cout << S << std::endl;
    
    
    S.propagate();
    
    std::cout << "\npropag:\n";
    std::cout << S << std::endl;
    
    S.restoreState(0);
    
    std::cout << "\nrestore:\n";
    std::cout << S << std::endl;
    
 
    
 
}
