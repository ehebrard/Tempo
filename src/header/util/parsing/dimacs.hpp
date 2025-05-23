#ifndef __DIMACSREADER_HH
#define __DIMACSREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Model.hpp"

namespace dimacs {

template <typename M, typename V> void parse(const std::string &fn, M &model, V& vars, bool& consistent) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");
    int n, l;
    size_t m;
    size_t clause_added{0};
    std::string dump;
//    std::vector<tempo::BooleanVar<int>> vars;
    std::vector<tempo::Literal<int>> cl;
    for (std::string line; getline(ifs, line);) {
      if (line[0] == 'c')
        continue;
      std::istringstream iss(line);
      if (line[0] == 'p') {
        iss >> dump;
        assert(dump == "p");
        iss >> dump;
        assert(dump == "cnf");
        iss >> n;
        iss >> m;
        for (auto i{0}; i < n; ++i) {
          auto x{model.newBoolean()};
          //          assert(x.id() == vars.size());
          vars.push_back(x);
        }
      } else {
        bool emptyline{false};
        while (true) {
          if (iss.eof()) {
            emptyline = true;
            break;
          }
          iss >> l;
          if (l == 0)
            break;
          if (l > 0) {
            cl.push_back(vars[l - 1] == true);
          } else {
            cl.push_back(vars[-l - 1] == false);
          }
        }
        if (not emptyline) {
            try {
                model.clauses.add(cl.begin(), cl.end());
                cl.clear();
                if (++clause_added == m)
                    break;
            } catch(tempo::Failure<int> &f) {
                consistent = false;
                return;
            }
        }
      }
    }
    //      assert(m == model.clauses.size());
  } catch (std::exception &e) {
    std::cout.flush();
    cerr << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
}

} // namespace dimacs

#endif
