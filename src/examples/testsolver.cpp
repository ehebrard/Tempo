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
#include <vector>


#include "Solver.hpp"


using namespace tempo;

void test1(Options& opt) {

    Solver<float> S(opt);
    
    std::cout << S << std::endl;
    
    
    auto b0{S.newBoolean()};
    auto b1{S.newBoolean()};
    auto b2{S.newBoolean()};

    
    auto x0{S.newTemporal()};
    auto x1{S.newNumeric()};
    auto x2{S.newTemporal()};
    auto x3{S.newTemporal()};
    
//    auto d0{S.newDisjunct(x0.before(x2,100), x0.after(x2,3))};
    auto d0{S.newDisjunct(x0.before(x2,100), x2.before(x0,3))};
    
    
    std::cout << S << std::endl;
    S.propagate();
    auto s1{S.saveState()};

    S.set(x0 <= float(18.999));
    
    S.set(x0 > float(-100));
    
    S.set(b1 == false);
    
    S.set(x1 < float(1000));
    
    S.set(x1 >= float(100));
    
    S.set(b2 == true);
    
    std::cout << S << std::endl;
    S.propagate();
    auto s2{S.saveState()};
    
    S.set(b0 == true);
    
    S.set(x1 <= float(300));
 
    S.set(x0 >= float(0));
    
    std::cout << S << std::endl;
    S.propagate();
    auto s3{S.saveState()};
    
    S.set(d0 == true);
    
    std::cout << S << std::endl;
    
    std::cout << "restore to state " << s3 << std::endl;
    S.restoreState(s3);
    std::cout << S << std::endl;
    
    S.propagate();
    auto s4{S.saveState()};
    
    S.set(d0 == false);
    
    std::cout << S << std::endl;
    
    std::cout << "restore to state " << s4 << std::endl;
    S.restoreState(s4);
    std::cout << S << std::endl;
    
    std::cout << "restore to state " << s2 << std::endl;
    S.restoreState(s2);
    std::cout << S << std::endl;
    
    std::cout << "restore to state " << s1 << std::endl;
    S.restoreState(s1);
    std::cout << S << std::endl;
    
}

void test2(Options& opt) {
    
    Solver<> S(opt);
    
    
    auto schedule{S.newJob(0,200)};
    
    S.set(schedule.start.after(0));
    S.set(schedule.start.before(0));
    
    
    Job<int> j0{S.newJob(15,15)};
    auto j1{S.newJob(7,10)};
    auto j2{S.newJob(12,12)};
    auto j3{S.newJob(1,5)};
    auto j4{S.newJob(10,10)};
    auto j5{S.newJob(6,6)};
    auto j6{S.newJob(10,12)};
    auto j7{S.newJob(7,7)};
    auto j8{S.newJob(3,3)};
    
    S.set(j0.start.after(schedule.start));
    S.set(j0.end.before(j1.start));
    S.set(j1.end.before(j2.start,10));
    S.set(j2.end.before(schedule.end));
    
    S.set(j3.start.after(schedule.start));
    S.set(j3.end.before(j4.start));
    S.set(j4.end.before(j5.start));
    S.set(j5.end.before(schedule.end));
    
    S.set(j6.start.after(schedule.start));
    S.set(j6.end.before(j7.start));
    S.set(j7.end.before(j8.start));
    S.set(j8.end.before(schedule.end));
    
    
    std::vector<Job<int>> jobsR1;
    jobsR1.emplace_back(j0);
    jobsR1.emplace_back(j3);
    jobsR1.emplace_back(j6);
    
    std::vector<Job<int>> jobsR2;
    jobsR2.emplace_back(j1);
    jobsR2.emplace_back(j4);
    jobsR2.emplace_back(j7);
    
    std::vector<Job<int>> jobsR3;
    jobsR3.emplace_back(j2);
    jobsR3.emplace_back(j5);
    jobsR3.emplace_back(j8);

    DisjunctiveResource<int> R1(S,jobsR1.begin(), jobsR1.end());
    DisjunctiveResource<int> R2(S,jobsR2.begin(), jobsR2.end());
    DisjunctiveResource<int> R3(S,jobsR3.begin(), jobsR3.end());

    std::vector<DisjunctVar<int>> X;
    
    R1.createOrderVariables(X);
    R2.createOrderVariables(X);
    R3.createOrderVariables(X);
    
    std::cout << S << std::endl;
    
}


int main(int argc, char *argv[]) {
    
    Options opt = tempo::parse(argc, argv);
    
    test2(opt);
}
