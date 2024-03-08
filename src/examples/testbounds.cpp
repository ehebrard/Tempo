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
    
    int n = 5;
    
    Options opt = tempo::parse(argc, argv);
    Scheduler<int> s_(opt);
    
    BoundSystem<int>& S{s_.domain.bounds};
    
    S.resize(n);
    
    std::cout << S << std::endl;
    
    for(auto l{0}; l<10; ++l) {
        
        ReversibleObject::env->save();
        
        try {
            for(auto i{0}; i<3; ++i) {
                
                event a{static_cast<event>(tempo::random()%n)};
                int k{static_cast<int>(tempo::random()%10000)};
                
                event b{static_cast<event>(tempo::random()%n)};
                int l{static_cast<int>(tempo::random()%100)};
                
                std::cout << "set " << prettyEvent(a) << " <= " << k << std::endl;
                
                S.set(UPPER, a, k);
                
                std::cout << "set " << prettyEvent(b) << " >= " << l << std::endl;
                
                S.set(LOWER, b, -l);
            }
            std::cout << std::endl << S << std::endl;
            
        } catch(std::exception& e) {
            std::cout << "negative cycle:\n" << S << std::endl;
            break;
        }
    }
    
    while(ReversibleObject::env->level() > 0) {
        ReversibleObject::env->restore(ReversibleObject::env->level()-1);
        std::cout << std::endl << S << std::endl;
    }
  
 
}
