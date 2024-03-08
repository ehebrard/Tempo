

#ifndef _TEMPO_CLAUSE_HPP
#define _TEMPO_CLAUSE_HPP

#include <vector>

#include "Global.hpp"
#include "Explanation.hpp"

namespace tempo {

class Clause : public std::vector<lit> { //}, public Explainer {

public:

  const static Clause* Empty;

  Clause(const int i = -1);
  virtual ~Clause();

  int id;
  index_t watcher_index[2];

  lit watcher(const bool) const;
  bool watch_rank(const lit) const;
size_t watch_index(const bool) const;

  bool operator<(const Clause &) const;
  bool operator==(const Clause &) const;
    
//    // explanation stuff
//     void xplain(const lit l, const hint h, std::vector<lit> &Cl) ;
//     std::ostream &print_reason(std::ostream &, const hint) const;
//     int getType() const;

  std::ostream& display(std::ostream &) const;

};

std::ostream &operator<<(std::ostream &, const Clause &);

}

#endif // _TEMPO_CLAUSE_HPP

