/************************************************
 * Tempo Profile.hpp
 * Implementation of the profile data structure as described in
 * Vincent Gingras and Claude-Guy Quimper. 2016. Generalizing the edge-finder rule for the cumulative constraint. In Proceedings of the Twenty-Fifth International Joint Conference on Artificial Intelligence (IJCAI'16). AAAI Press, 3103–3109.
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/



#ifndef _TEMPO_LIST_HPP
#define _TEMPO_LIST_HPP

#include <vector>

namespace tempo {

template<typename E>
struct Node {

  Node();

  template <typename... T>
  Node(const int n, const int p, T... args)
      : next(n), prev(p), content(args...) {}

  int next{0};
  int prev{0};
  E content;
};





template<typename E>
class List {
    
    
     class iterator
    {
        friend class List<E>;
        
    public:
        iterator(List<E>* _l_, const int _i_) : l(_l_), index(_i_) {}
        
        const E& operator*() const {
            return l->nodes[index].content;
        }
        
        void operator++() {
            index = l->nodes[index].next;
        }
        
        void operator--() {
            index = l->nodes[index].prev;
        }
        
        bool operator==(const List<E>::iterator li) {
            return li.l == l and li.index == index;
        }
        
        bool operator!=(const List<E>::iterator li) {
            return not this->operator==(li);
        }
        
        E* operator->() {
               return &(l->nodes[index].content);
           }

        int next() { return l->nodes[index].next; }

        List<E>* l;
        int index;
 
    };
    
    
    class reverse_iterator
   {
       friend class List<E>;
       
   public:
       reverse_iterator(List<E>* _l_, const int _i_) : l(_l_), index(_i_) {}
       
       const E& operator*() const {
           return l->nodes[index].content;
       }
       
       void operator++() {
           index = l->nodes[index].prev;
       }
       
       void operator--() {
           index = l->nodes[index].next;
       }
       
       bool operator==(const List<E>::reverse_iterator li) {
           return li.l == l and li.index == index;
       }
       
       bool operator!=(const List<E>::reverse_iterator li) {
           return not this->operator==(li);
       }
       
       E* operator->() {
              return &(l->nodes[index].content);
          }

       int next() { return l->nodes[index].prev; }

       List<E>* l;
       int index;

   };

public:
  List() {
    nodes.resize(1); // dummy element to represent head/tail
  }

    bool empty() const {return nodes[tail].next == tail;}

    void pop_back() {
      remove(nodes.size() - 1);
      nodes.pop_back();
    }

    template <typename... T>
    int create_element(T ...args) {
        int tindex{static_cast<int>(nodes.size())};
        nodes.emplace_back(-1,-1,args...);
        return tindex;
    }

    int first() const { return nodes[tail].next; }
    int last() const { return nodes[tail].prev; }

    //    int tail() const {return nodes[tail].prev;}

    int next(const int i) const { return nodes[i].next; }
    int prev(const int i) const { return nodes[i].prev; }

    void add_front(const int tindex) { add_after(tail, tindex); }

    void add_after(const int rid, const int tindex) {
        nodes[tindex].prev = rid;
        int nid = nodes[rid].next;
        nodes[rid].next = tindex;
        nodes[nid].prev = tindex;
        nodes[tindex].next = nid;
    }

    void add_before(const int rid, const int tindex) {
      nodes[tindex].next = rid;
      int pid = nodes[rid].prev;
      nodes[rid].prev = tindex;
      nodes[pid].next = tindex;
      nodes[tindex].prev = pid;
    }

    template <typename Func> void add_when(Func test, const int tindex) {
      for (auto i{begin()}; i != end(); ++i) {
        if (test(i))
          add_after(i.index, tindex);
      }
    }

    void remove(const int tindex) {
      nodes[nodes[tindex].prev].next = nodes[tindex].next;
      nodes[nodes[tindex].next].prev = nodes[tindex].prev;
    }

    //    void cut_interval(const int i, const int j) {

    void clear() {
      nodes[tail].next = tail;
      nodes[tail].prev = tail;
    }

    List<E>::iterator begin() { return iterator(this, first()); }
    List<E>::iterator end() { return iterator(this, tail); }
    
    List<E>::reverse_iterator rbegin() { return reverse_iterator(this, last()); }
    List<E>::reverse_iterator rend() { return reverse_iterator(this, tail); }
    
    
//    List<E>::const_iterator begin() const { return const_iterator(this, nodes[tail].next); }
//    List<E>::const_iterator end() const { return const_iterator(this, tail); }
    
    List<E>::iterator at(const int i) { return iterator(this, i); }
    
    size_t size() const {return nodes.size()-1;}
    E& operator[](const int i) {return nodes[i].content;}

    std::ostream &display(std::ostream &os) {
      for (auto tp{begin()}; tp != end(); ++tp) {
        os << tp.index << ": " << *tp << std::endl;
      }
      return os;
    }

    static int tail;

  protected:
    std::vector<Node<E>> nodes;
};

template<typename E>
int List<E>::tail = 0;

template <typename E> Node<E>::Node() {
  next = List<E>::tail;
  prev = List<E>::tail;
}

template<typename E>
std::ostream &operator<<(std::ostream &os, List<E> &x) {
    return x.display(os);
}

} // namespace tempo

#endif // _TEMPO_THETATREE_HPP
