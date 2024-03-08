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


#include "Global.hpp"
#include "DirectedGraph.hpp"


using namespace tempo;




int main() {
    
    DirectedGraph<int> g(10);
    
    
        for(auto l{0}; l<5; ++l) {
    
            ReversibleObject::env->save();
    
            for(auto l{0}; l<3; ++l) {
                auto a{tempo::random() % 10};
                auto b{tempo::random() % 10};
                if(a != b)
                {
                    g.add(a,b);
                }
            }
    
            std::cout << g << std::endl;
        }

    
    for(auto l{0}; l<3; ++l) {
        
        
        auto a{tempo::random() % 10};
        auto b{tempo::random() % 10};
 
        if(a != b and g.has(a) and g.has(b)) {
            
            ReversibleObject::env->save();
            
            std::cout << "merge " << b << " to " << a << std::endl;
            
            g.merge(a,b);
            
            std::cout << g << std::endl;
            
        }
        
    }
    
    for(auto l{0}; l<3; ++l) {
        
        
        auto a{tempo::random() % 10};
 
        if(g.has(a)) {
            
            ReversibleObject::env->save();
            
            std::cout << "remove " << a << std::endl;
            
            g.remove(a);
            
            std::cout << g << std::endl;
            
        }
        
    }
     
    while(ReversibleObject::env->level() > 0) {
        ReversibleObject::env->restore(ReversibleObject::env->level()-1);
        
        std::cout << g << std::endl;
    }
    
    
    
    
    
    
    DirectedGraph<StampedLabeledEdge<float>> gl(10);

    for(auto l{0}; l<5; ++l) {
        
        ReversibleObject::env->save();
        
        for(auto l{0}; l<3; ++l) {
            auto a{tempo::random() % 10};
            auto b{tempo::random() % 10};
            float w{static_cast<float>(tempo::random() % 100)/static_cast<float>(1+tempo::random() % 10)};
            if(a != b)
            {
                gl.add(a,StampedLabeledEdge<float>(b,w));
            }
        }
        
        std::cout << gl << std::endl;
    }
    
    for(auto l{0}; l<3; ++l) {
        
        
        auto a{tempo::random() % 10};
        auto b{tempo::random() % 10};
 
        if(a != b and gl.has(a) and gl.has(b)) {
            
            ReversibleObject::env->save();
            
            std::cout << "merge " << b << " to " << a << std::endl;
            
            gl.merge(a,b);
            
            std::cout << gl << std::endl;
            
        }
        
    }
    
    for(auto l{0}; l<3; ++l) {
        
        
        auto a{tempo::random() % 10};
 
        if(gl.has(a)) {
            
            ReversibleObject::env->save();
            
            std::cout << "remove " << a << std::endl;
            
            gl.remove(a);
            
            std::cout << gl << std::endl;
            
        }
        
    }
     
    while(ReversibleObject::env->level() > 0) {
        ReversibleObject::env->restore(ReversibleObject::env->level()-1);
        
        std::cout << gl << std::endl;
    }
    
    
    
    
    
 
    
}
