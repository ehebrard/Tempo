#ifndef __DIMACSREADER_HH
#define __DIMACSREADER_HH

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Model.hpp"

namespace dimacs {

template <typename M, typename V> void parse(const std::string &fn, M &model, V& vars) {
  using std::cerr;
  try {
    std::ifstream ifs(fn);
    if (!ifs)
      throw std::runtime_error("Could not open file for reading");
    int n, l;
    size_t m;
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
          assert(x.id() == vars.size());
          vars.push_back(x);
        }
      } else {
        while (true) {
          iss >> l;
          if (l == 0)
            break;
          if (l > 0) {
            cl.push_back(vars[l - 1] == true);
          } else {
            cl.push_back(vars[-l - 1] == false);
          }
        }
        model.clauses.add(cl.begin(), cl.end());
        cl.clear();
        if (model.clauses.size() == m)
          break;
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
