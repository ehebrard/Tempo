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


#include "TemporalNetwork.hpp"


using namespace tempo;




int main() {
    
    size_t n{10};
    
    TemporalNetwork<int> N;
    
    N.resize(n);
    
    while(true) {
        
        ReversibleObject::env->save();
        
        event x{static_cast<event>(tempo::random()%n)};
        event y{static_cast<event>(tempo::random()%n)};
        
        int d{static_cast<int>(tempo::random()%100)};
        if(y == 0) {
            d = -d;
        } else {
            if(x != 1) {
                if((tempo::random()%2)==1)
                    d = -d;
            } else {
                if(y == 0)
                    d += 300;
            }
        }
          
        try {
            N.newEdge(x,y,d);
        } catch(std::exception &e) {
            std::cout << e.what() << std::endl;
            std::cout << N << std::endl;
            break;
        }
        std::cout << N << std::endl;
    }
    
    
    
    
    while(ReversibleObject::env->level() > 0) {
        ReversibleObject::env->restore(ReversibleObject::env->level()-1);
        
        std::cout << N << std::endl;
    }
    
 
}
