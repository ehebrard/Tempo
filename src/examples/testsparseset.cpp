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

#include "SparseSet.hpp"


using namespace tempo;





int main() {
    
    SparseSet<long,size_t> s(100);
    
    std::cout << s << std::endl;
    
    
    for(auto i{0}; i<50; ++i) {
        auto r{random() % 100};
        if(not s.has(r))
            s.add(r);
    }
    
    std::cout << s << std::endl;
    
    s.fill();
    
    std::cout << s << std::endl;
    
    s.clear();
    
    std::cout << s << std::endl;
    
    
}
