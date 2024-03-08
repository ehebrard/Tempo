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
for more detailrevset.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/


#include <iostream>

#include "ReversibleObject.hpp"
#include "SparseSet.hpp"


using namespace tempo;



//int main(int argc, char *argv[]) {
int main() {
    
    
    
    
//    ReversibleVector<int> revvec({1,7,12,34});
////    ReVector<int> revvec({1,7,12,34});
//
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->save();
//    std::cout << "SAVE\n";
//    
//    revvec.push_back(0);
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->save();
//    std::cout << "SAVE\n";
//    
//    revvec.push_back(1);
//    revvec.push_back(2);
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->save();
//    std::cout << "SAVE\n";
//    
//    revvec.push_back(3);
//    revvec.push_back(4);
//    revvec.push_back(5);
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->restore(ReversibleObject::env->level()-1);
//    std::cout << "RESTORE\n";
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->restore(ReversibleObject::env->level()-1);
//    std::cout << "RESTORE\n";
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
//    
//    ReversibleObject::env->restore(ReversibleObject::env->level()-1);
//    std::cout << "RESTORE\n";
//    
//    for(auto x : revvec)
//        std::cout << " " << x;
//    std::cout << std::endl;
    
    
    Reversible<int> revint{0};
    SparseSet<Reversible<size_t>> revset(20);
    
    BacktrackEnvironment * another_env = new BacktrackEnvironment();
    ReversibleVector<int> revvec({1,7,12,34}, another_env);

    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    ReversibleObject::env->save();
    std::cout << "SAVE\n";
    
    revint = 13;
    for(auto i{0}; i<5; ++i) {
        auto r{tempo::random() % 20};
        if(not revset.has(r))
            revset.add(r);
    }
    revvec.push_back(0);
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    another_env->save();
    std::cout << "SAVE*\n";
    
    revint = 7;
    for(auto i{0}; i<5; ++i) {
        auto r{tempo::random() % 20};
        if(not revset.has(r)) {
            revset.add(r);
        }
    }
    revvec.push_back(1);
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    ReversibleObject::env->save();
    std::cout << "SAVE\n";
    
    revint = 9;
    for(auto i{0}; i<5; ++i) {
        auto r{tempo::random() % 20};
        if(not revset.has(r)) {
            revset.add(r);
        }
    }
    revvec.push_back(2);
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    another_env->save();
    std::cout << "SAVE*\n";
    
    revint = 28;
    for(auto i{0}; i<5; ++i) {
        auto r{tempo::random() % 20};
        if(not revset.has(r)) {
            revset.add(r);
        }
    }
    revvec.push_back(3);
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    
    ReversibleObject::env->restore(ReversibleObject::env->level()-1);
    std::cout << "RESTORE\n";
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    
    another_env->restore(another_env->level()-1);
    std::cout << "RESTORE*\n";
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    
    ReversibleObject::env->restore(ReversibleObject::env->level()-1);
    std::cout << "RESTORE\n";
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    
    another_env->restore(another_env->level()-1);
    std::cout << "RESTORE*\n";
    
    std::cout << revint << std::endl;
    std::cout << revset << std::endl;
    for(auto x : revvec)
        std::cout << " " << x;
    std::cout << std::endl;
    
    
}
