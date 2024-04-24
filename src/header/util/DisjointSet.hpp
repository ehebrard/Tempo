
#ifndef _TEMPO_DISJOINTSET_HPP
#define _TEMPO_DISJOINTSET_HPP

#include <iostream>
#include <vector>

#include "ReversibleObject.hpp"

namespace tempo {


/**********************************************
 * DisjointSet
 **********************************************/
/// Disjoint set representation


template<typename E=int>
class DisjointSet {
   
private:
    
  std::vector<E> parent;
  std::vector<uin16_t> size;

    
public:
  /*!@name Constructors*/
  //@{
  explicit DisjointSet(const size_t n = 0);
  void resize(const size_t n);

  /*!@name Accessors*/
  //@{
  size_t find(const E elt);
    void union(const E x, const E y);
    
    // to use merge, x and y must be roots !!
    void merge(const E x, const E y);
  //@}


  /*!@name Miscellaneous*/
  //@{
  std::ostream &display(std::ostream &os) const;
    //@}
};


template<typename E>
DisjointSet<E>::DisjointSet(const size_t n) {
    resize(n);
}


template<typename E>
void DisjointSet<E>::resize(const size_t n) {
  while (parent.size() < n) {
      parent.push_back(static_cast<E>(parent.size()));
  }
    size.resize(n,1);
}

template<typename E>
size_t DisjointSet<E>::find(const E x) {
    while(parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x
}

template<typename E>
void DisjointSet<E>::union(const E x, const E y) {
    auto sx{find(x)};
    auto sy{find(y)};
    if(sx != sy) {
        merge(sx, sy);
    }
}

template<typename E>
void DisjointSet<E>::merge(const E x, const E y) {
    if(size[x] > size[y]) {
        parent[y] = x;
        size[x] += size[y];
    } else {
        parent[x] = y;
        size[y] += size[x];
    }
}


template<typename E>
std::ostream &DisjointSet<E>::display(std::ostream &os) const {
    for(size_t i{0}; i<parent.size(); ++i) {
        os << setw(4) << i;
    }
    os << std::endl;
    for(auto p : parent) {
        os << setw(4) << p;
    }
    os << std::endl;
    for(auto s : size) {
        os << setw(4) << s;
    }
    os << std::endl;
}

template<typename E>
std::ostream &operator<<(std::ostream &os, const DisjointSet<E,T> &x) {
  return x.display(os);
}

template<typename E>
std::ostream &operator<<(std::ostream &os, const DisjointSet<E,T> *x) {
  return (x ? x->display(os) : os);
}

}


#endif // _TEMPO_DisjointSet_HPP
