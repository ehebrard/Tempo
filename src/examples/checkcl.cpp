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


#include <chrono>
#include <iostream>
#include <vector>
#include <filesystem>
#include <optional>

#include "Scheduler.hpp"
#include "util/parsing/format.hpp"
#include "util/parsing/jsp.hpp"
#include "util/parsing/osp.hpp"


using namespace tempo;


int main(int argc, char *argv[]) {
  Options opt = tempo::parse(argc, argv);

  seed(opt.seed);

  ProblemInstance data;


  if (opt.input_format == "osp") {
    data = osp::read_instance(opt.instance_file);
  }
  else if (opt.input_format == "jsp") {
    data = jsp::read_instance(opt.instance_file);
  }

  Scheduler<int> S(opt);
  for (auto d : data.durations) {
    S.newTask(d, d);
  }
    
  for (auto [x, y, d] : data.constraints) {

      assert(d==0);

    S.newPrecedence(x, y, d);

  }
    
  for (auto &job : data.resources) {
      for (auto ti{job.begin()}; ti!=job.end(); ++ti) {
          for (auto tj{ti+1}; tj!=job.end(); ++tj) {
              S.newVariable({START(*ti), END(*tj), 0}, {START(*tj), END(*ti), 0});
          }
      }
  }
    
//    S.display(std::cout, true, true);

//    std::cout << "\nsave:\n" << S << std::endl;
    
    std::ifstream cl_file("cl.txt", std::ifstream::in);
    
    int t, n, x, y, d;
    
    std::vector<int> X;
    std::vector<int> Y;
    std::vector<int> D;
    
    int line{0};
    do {
        
        X.clear();
        Y.clear();
        D.clear();
        
        S.saveState();
        
        cl_file >> t;
        cl_file >> n;
        if(not cl_file.good())
            break;
        
        bool trivially_unsat{false};
        
        
        
        
        for(auto i{0}; i<n; ++i) {
            cl_file >> x;
            cl_file >> y;
            cl_file >> d;
            
            X.push_back(x);
            Y.push_back(y);
            D.push_back(d);
            
//            std::cout << "add " << prettyEvent(y) << " - " << prettyEvent(x) << " <= " << d << std::endl;
            
            try{
                S.newMaximumLag(x, y, d);
            } catch(Failure& f) {
                trivially_unsat = true;
            }
        }
        
//        std::cout << "\nsolve:\n" << S << std::endl;
        
        
        if(not trivially_unsat)
        {
//            std::cout << "ok (trivial)\n";
//        } else {
            
//            std::cout << S << std::endl;
//            exit(1);
            
            S.search();
            
            if(S.satisfiable()) {
                std::cout << "cl " << line << " (" << (t ? "expl" : "cut") << "): ";
                std::cout << "bug!\n";
                
                for(size_t i{0}; i<X.size(); ++i) {
                    std::cout << "> " << prettyEvent(Y[i]) << " - " << prettyEvent(X[i]) << " <= " << D[i] << std::endl;
                }
                
                exit(1);
            } 
//            else {
//                std::cout << "ok\n";
//            }
        }
        ++line;
        
//        std::cout << S.env.level() << " --> 0"
        
        S.restoreState(0);
        
//        std::cout << "\nrestore:\n" << S << std::endl;
        
    } while(true);
    
    
    
 
 
}
