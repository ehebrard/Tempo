
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
  std::vector<uint16_t> size;

    
public:
  /*!@name Constructors*/
  //@{
  explicit DisjointSet(const size_t n = 0);
  void resize(const size_t n);

  /*!@name Accessors*/
  //@{
  size_t find(const E elt);
    void merge(const E x, const E y);
    void clear();
    template<typename Iter>
    void init(Iter, Iter);
    // to use merge_roots, x and y must be roots !!
    void merge_roots(const E x, const E y);
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
void DisjointSet<E>::clear() {
    for(size_t i{0}; i<parent.size(); ++i) {
        parent[i] = i;
        size[i] = 1;
    }
}

template<typename E>
template<typename Iter>
void DisjointSet<E>::init(Iter beg_elt, Iter end_elt) {
    for(auto it{beg_elt}; it!=end_elt; ++it) {
        parent[*it] = *it;
        size[*it] = 1;
    }
}

template<typename E>
size_t DisjointSet<E>::find(const E elt) {
    auto x{elt};
    while(parent[x] != x) {
        
        std::cout << "p[" << x << "] <- p[" << parent[x] << "]\n";
        
        
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

template<typename E>
void DisjointSet<E>::merge(const E x, const E y) {
    auto sx{find(x)};
    auto sy{find(y)};
    if(sx != sy) {
        merge_roots(sx, sy);
    }
}

template<typename E>
void DisjointSet<E>::merge_roots(const E x, const E y) {
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
        os << std::setw(4) << i;
    }
    os << std::endl;
    for(auto p : parent) {
        os << std::setw(4) << p;
    }
    os << std::endl;
    for(auto s : size) {
        os << std::setw(4) << s;
    }
    os << std::endl;
    return os;
}

template<typename E>
std::ostream &operator<<(std::ostream &os, const DisjointSet<E> &x) {
  return x.display(os);
}

template<typename E>
std::ostream &operator<<(std::ostream &os, const DisjointSet<E> *x) {
  return (x ? x->display(os) : os);
}

}


#endif // _TEMPO_DisjointSet_HPP
