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
    
    Node() {}
    
    template <typename... T>
    Node(const int n, const int p, T ...args) : next(n), prev(p), content(args...) {}
    
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
        
        List<E>* l;
        int index;
 
    };

public:
    
    List() {nodes.resize(1);}
    
    bool empty() const {return nodes[tail].next == tail;}
    
    
    template <typename... T>
    int create_element(T ...args) {
        int tindex{static_cast<int>(nodes.size())};
        nodes.emplace_back(-1,-1,args...);
        return tindex;
    }
    
    int next(const int i) const {return nodes[i].next;}
    
    void add_front(const int tindex) {
//        nodes[tindex].prev = tail;
//        nodes[tindex].next = nodes[tail].next;
//        nodes[nodes[tail].next].prev = tindex;
//        nodes[tail].next = tindex;
        add_after(tail, tindex);
    }
    
    void add_after(const int rid, const int tindex) {
        nodes[tindex].prev = rid;
        int nid = nodes[rid].next;
        nodes[rid].next = tindex;
        nodes[nid].prev = tindex;
        nodes[tindex].next = nid;
    }
    
    void remove(const int tindex) {
//        
//        std::cout << "rm " << tindex << std::endl;
//        std::cout << "next = " << nodes[tindex].next << std::endl;
//        std::cout << "prev = " << nodes[tindex].prev << std::endl;
//        
//        
        nodes[nodes[tindex].prev].next = nodes[tindex].next;
        nodes[nodes[tindex].next].prev = nodes[tindex].prev;
    }
    
    List<E>::iterator begin() { return iterator(this, nodes[tail].next); }
    List<E>::iterator end() { return iterator(this, tail); }
    
    
//    List<E>::const_iterator begin() const { return const_iterator(this, nodes[tail].next); }
//    List<E>::const_iterator end() const { return const_iterator(this, tail); }
    
    List<E>::iterator at(const int i) { return iterator(this, i); }
    
    size_t size() const {return nodes.size()-1;}
    E& operator[](const int i) {return nodes[i].content;}

    std::ostream &display(std::ostream &os) {
        

        
//        os << "head = " << next[tail];
//        os << "\nnext:";
//                for(auto n : next) {
//                    os << " " << n;
//                }
//                os << std::endl << "prev:";
//                for(auto p : prev) {
//                    os << " " << p;
//                }
//        os << std::endl;
        
        for(auto tp{begin()}; tp!=end(); ++tp) {
            os << tp.index << ": " << *tp << std::endl;
        }
        return os;
    }
    
protected:
    
    std::vector<Node<E>> nodes;
    
    static int tail;
};



template<typename E>
int List<E>::tail = 0;


template<typename E>
std::ostream &operator<<(std::ostream &os, List<E> &x) {
    return x.display(os);
}



//
//template<typename E>
//class Machin {
//    
//    
//     class iterator
//    {
//        friend class Machin<E>;
//        
//    public:
//        iterator(const Machin<E>* _l_, const int _i_) : l(_l_), index(_i_) {}
//        
//        const E& operator*() const {
//            return l->elements[index];
//        }
//        
//        void operator++() {
//            index = l->next[index];
//        }
//        
//        void operator--() {
//            index = l->prev[index];
//        }
//        
//        bool operator==(const Machin<E>::iterator li) {
//            return li.l == l and li.index == index;
//        }
//        
//        bool operator!=(const Machin<E>::iterator li) {
//            return not this->operator==(li);
//        }
//        
//        const Machin<E>* l;
//        int index;
// 
//    };
//
//public:
//    
//    Machin() {next.resize(1,0); prev.resize(1,0); elements.resize(1);}
//    
//    bool empty() const {return next[tail] == tail;}
//    
//    
//    template <typename... T>
//    int create_element(T ...args) {
//        int tindex{static_cast<int>(elements.size())};
//        elements.emplace_back(args...);
//        
//        prev.resize(elements.size(), -1);
//        next.resize(elements.size(), -1);
//        
//        return tindex;
//    }
//    
//    void add_front(const int tindex) {
//        next[tindex] = next[tail];
//        prev[next[tail]] = tindex;
//        next[tail] = tindex;
//    }
//    
//    void add_after(const int rid, const int tindex) {
//        prev[tindex] = rid;
//        int nid = next[rid];
//        next[rid] = tindex;
//        
//        prev[nid] = tindex;
//        next[tindex] = nid;
//    }
//    
//    void remove(const int tindex) {
//        next[prev[tindex]] = next[tindex];
//        prev[next[tindex]] = prev[tindex];
//    }
//    
//    Machin<E>::iterator begin() const { return iterator(this, next[tail]); }
//    Machin<E>::iterator end() const { return iterator(this, tail); }
//    
//    Machin<E>::iterator at(const int i) const { return iterator(this, i); }
//    
//    size_t size() const {return elements.size()-1;}
//    E& operator[](const int i) {return elements[i];}
//
//    std::ostream &display(std::ostream &os) const {
//        
//
//        
////        os << "head = " << next[tail];
////        os << "\nnext:";
////                for(auto n : next) {
////                    os << " " << n;
////                }
////                os << std::endl << "prev:";
////                for(auto p : prev) {
////                    os << " " << p;
////                }
////        os << std::endl;
//        
//        for(auto tp{begin()}; tp!=end(); ++tp) {
//            os << tp.index << ": " << *tp << std::endl;
//        }
//        return os;
//    }
//    
//protected:
//    
//    std::vector<E> elements;
//    std::vector<int> next;
//    std::vector<int> prev;
//    
//    static int tail;
// 
////    int head{0};
//};
//
//template<typename E>
//int Machin<E>::tail = 0;
//
//
//template<typename E>
//std::ostream &operator<<(std::ostream &os, const Machin<E> &x) {
//    return x.display(os);
//}

} // namespace tempo

#endif // _TEMPO_THETATREE_HPP
