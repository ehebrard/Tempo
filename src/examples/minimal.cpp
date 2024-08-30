#include <iostream>
#include <vector>

#include "Solver.hpp"

using namespace tempo;


int main() {

  tempo::Options o = tempo::no_option;
  Solver<> solver(o);
 
    auto horizon{solver.newNumeric()};
    
    auto schedule{solver.between(solver.zero(), horizon)};
    
    auto s1{solver.newNumeric()};
    auto t1{solver.between(s1, s1+20)};
    
    
    auto s2{solver.newNumeric()};
    auto t2{solver.between(s2, s2+10)};
    
    
    auto s3{solver.newNumeric()};
    auto t3{solver.between(s3, s3+30)};
    
    
    auto s4{solver.newNumeric()};
    auto t4{solver.between(s4, s4+25)};
    
    
    auto s5{solver.newNumeric()};
    auto t5{solver.between(s5, s5+15)};
    
    
    solver.post(s1 >= 0);
    solver.post(s2 >= 0);
    solver.post(s3 >= 0);
    solver.post(s4 >= 0);
    solver.post(s5 >= 0);
    
    solver.post(t1.end <= schedule.end);
    solver.post(t2.end <= schedule.end);
    solver.post(t3.end <= schedule.end);
    solver.post(t4.end <= schedule.end);
    solver.post(t5.end <= schedule.end);
    
    
    solver.post(NoOverlap<int>(schedule, {t1,t2,t3}));

    solver.post(NoOverlap<int>(schedule, {t4, t5}));
    
    std::cout << solver << std::endl;
    
    solver.minimize(schedule.duration);
    
    

}
