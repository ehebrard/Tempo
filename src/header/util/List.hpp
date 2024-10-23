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
      last_saved = tail;
      
#ifdef DBG_LIST
      verify("constructor");
#endif
  }

    bool empty() const {return nodes[tail].next == tail;}

    void pop_back() {
      remove_and_forget(nodes.size() - 1);
      nodes.pop_back();
//        present.pop_back();
        
#ifdef DBG_LIST
      verify("pop back");
#endif
    }

    template <typename... T>
    int create_element(T ...args) {
        int tindex{static_cast<int>(nodes.size())};
        nodes.emplace_back(-1,-1,args...);
//        present.push_back(false);
        
#ifdef DBG_LIST
      verify("create element");
#endif
        
        return tindex;
    }

    int first() const { return nodes[tail].next; }
    int last() const { return nodes[tail].prev; }

    int next(const int i) const { return nodes[i].next; }
    int prev(const int i) const { return nodes[i].prev; }

    void add_front(const int tindex) { add_after(tail, tindex);
#ifdef DBG_LIST
      verify("add front");
#endif
    }

    void add_after(const int rid, const int tindex) {
#ifdef DBG_LIST
        std::stringstream msg;
        msg << "before add " << tindex << " after " << rid ;
      verify(msg.str());
#endif
        
        ++add_count;
        
//        if(add_count == 3021138) {
//            std::cout << std::endl;
//            display(std::cout);
//        }
        
        
        nodes[tindex].prev = rid;
        int nid = nodes[rid].next;
        nodes[rid].next = tindex;
        nodes[nid].prev = tindex;
        nodes[tindex].next = nid;
        
//        if(nodes[tindex].next == tindex) {
//            std::cout << "cycle!! " << add_count << " " << tindex << " " << rid << "\n";
//            exit(1);
//        }
        
//        present[tindex] = true;
#ifdef DBG_LIST
        msg.clear();
        msg << "after add " << tindex << " after " << rid ;
      verify(msg.str());
#endif
    }

    void add_before(const int rid, const int tindex) {
      nodes[tindex].next = rid;
      int pid = nodes[rid].prev;
      nodes[rid].prev = tindex;
      nodes[pid].next = tindex;
      nodes[tindex].prev = pid;
//        present[tindex] = true;
        
#ifdef DBG_LIST
      verify("add before");
#endif
    }

    template <typename Func> void add_when(Func test, const int tindex) {
      for (auto i{begin()}; i != end(); ++i) {
        if (test(i))
          add_after(i.index, tindex);
      }
        
#ifdef DBG_LIST
      verify("add when");
#endif
    }

    void remove(const int tindex) {
        remove_and_forget(tindex);
        nodes[tindex].next = last_saved;
        last_saved = tindex;

#ifdef DBG_LIST
      verify("remove");
#endif
    }
    
    void remove_and_forget(const int tindex) {

#ifdef DBG_LIST
      verify("before remove and forget");
#endif

      nodes[nodes[tindex].prev].next = nodes[tindex].next;
      nodes[nodes[tindex].next].prev = nodes[tindex].prev;

#ifdef DBG_LIST
      verify("after remove and forget");
#endif
    }
    
    void re_add() {
        auto x{last_saved};
        last_saved = nodes[x].next;
        add_after(nodes[x].prev, x);
        
#ifdef DBG_LIST
      verify("re add");
#endif
    }
 
    void clear() {
      nodes[tail].next = tail;
      nodes[tail].prev = tail;
        last_saved = tail;
        
#ifdef DBG_LIST
      verify("clear");
#endif
    }

    List<E>::iterator begin() { return iterator(this, first()); }
    List<E>::iterator end() { return iterator(this, tail); }

    List<E>::reverse_iterator rbegin() { return reverse_iterator(this, last()); }
    List<E>::reverse_iterator rend() { return reverse_iterator(this, tail); }
    
    
//    List<E>::const_iterator begin() const { return const_iterator(this, nodes[tail].next); }
//    List<E>::const_iterator end() const { return const_iterator(this, tail); }
    
    List<E>::iterator at(const int i) { return iterator(this, i); }
    List<E>::iterator reverse_at(const int i) {
      return reverse_iterator(this, i);
    }

    size_t size() const {return nodes.size()-1;}
    E& operator[](const int i) {return nodes[i].content;}
    const E &operator[](const int i) const { return nodes[i].content; }

    std::ostream &display(std::ostream &os) {
#ifdef DBG_LIST
        int itermax{1000};
#endif
        
        auto width{static_cast<int>(std::log(nodes.size()))+1};
      for (auto tp{begin()}; tp != end(); ++tp) {
        os << std::setw(width) << tp.index << ": " << *tp << " " << prev(tp.index) << "|" << next(tp.index) << std::endl;
          
#ifdef DBG_LIST
        if(--itermax == 0)
            exit(1);
#endif
      }
      return os;
    }

#ifdef DBG_LIST
    void verify(std::string msg) {

        ++op_count;

        int i, p;
        
//        int itermax{1000};

        if(op_count >= 15440720) {
            std::cout << op_count << ": " << msg << std::endl;
            display(std::cout);
        }

        
        i = p = tail;
        do {
            i = next(i);
            if(prev(i) != p) {
                std::cout << msg << " (" << p << "->" << i << "<-" << prev(i) << ") @" << op_count << "\n" ;
//                display(std::cout);
                exit(1);
            }
            p = i;
            
//            if(--itermax == 0)
//            {
//                std::cout << "cycle?\n";
//                exit(1);
//            }
        } while(i != tail);
        
//        itermax = 1000;
        i = p = tail;
        do {
             i = prev(i);
            if(next(i) != p) {
                std::cout << msg << " (" << p << "<-" << i << "->" << next(i) << ") @" << op_count << "\n" ;
//                display(std::cout);
                exit(1);
            }
            p = i;
//            if(--itermax == 0)
//            {
//                std::cout << "cycle?\n";
//                exit(1);
//            }
        } while(i != tail);
    }
#endif
        

    static int tail;

  protected:
    std::vector<Node<E>> nodes;
//    std::vector<bool> present;
    int last_saved;
    
#ifdef DBG_LIST
    unsigned long op_count{0};
#endif
    
public:
    unsigned long add_count{0};
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
