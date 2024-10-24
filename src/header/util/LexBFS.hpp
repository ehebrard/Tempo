
#ifndef _TEMPO_LEXBFS_HPP
#define _TEMPO_LEXBFS_HPP

#include <iostream>
#include <vector>
#include <assert.h>

#include "Global.hpp"

//namespace tempo {


/**********************************************
 * LexBFS
 **********************************************/

using pos_t = tempo::index_t;
using vertex_t = tempo::index_t;
using interval_t = tempo::index_t;


struct VertexInfo {
    pos_t position; // in the global ordering vector
    interval_t interval;
    
    VertexInfo() {}
    VertexInfo(const pos_t p, const interval_t i) : position(p), interval(i) {}
};


struct IntervalInfo {
    pos_t first; // position of the first element in the global ordering vector
    vertex_t stamp; // pivot responsible for the creation of this interval
    
    IntervalInfo() {}
    IntervalInfo(const pos_t f, const vertex_t s) : first(f), stamp(s) {}
};


class LexBFS {
    
public:
    
    // contains all vertices, partitioned into sets
    std::vector<vertex_t> ordering;
    
    /*!@name Constructors*/
    //@{
    LexBFS();
    template<typename Set>
    void initialize(const size_t n, Set &S);
    //@}
    
    /*!@name Accessor*/
    //@{
    template<typename Set, typename Graph>
    void explore(Set &S, Graph &G);
    // whether x comes before y in the simplicial order (after in ordering)
    bool before(vertex_t x, vertex_t y) const;
    //@}
    
    /*!@name Helpers*/
    //@{
    // Find and remove a vertex v from the first interval. Remove the interval if it is empty
    void pop_front(interval_t interval);
    void remove(interval_t interval);
    void create_interval(interval_t interval);
    void clear();
    void move(vertex_t vertex);
    //@}
    
    
    /*!@name Miscellaneous*/
    //@{
    std::ostream &display(std::ostream &os) const;
    //@}
    
private:
    
    // contains the necessary information for each vertex: its position in 'vertices' and the set it belongs to
    std::vector<VertexInfo> vertices_info;
    
    // set of pointers to the start of an interval in 'vertices'
    std::vector<IntervalInfo> intervals_info;
    
    // intervals are ordered in a list
    std::vector<interval_t> next_interval;
    std::vector<interval_t> prev_interval;
    
    // the first interval
    static const interval_t head;
    // the last (dummy) interval
    static const interval_t tail;
    
    static const interval_t no_interval;
    
    static const interval_t default_interval;
    
    static const vertex_t no_vertex;
    
    // interval indices that have been released
    std::vector<interval_t> free_indices;
    
};



std::ostream &operator<<(std::ostream &os, const LexBFS &x);

template<typename Set>
void LexBFS::initialize(const size_t n, Set& V) {
    ordering.clear();
//    while(vertices_info.size() < n) {
//        int i{static_cast<int>(vertices_info.size())};
//        vertices_info.emplace_back(i,default_interval);
//    }
    vertices_info.resize(n, {0, no_interval});
    pos_t p{0};
    for(auto v : V) {
        ordering.push_back(v);
        vertices_info[v] = {p++, default_interval};
    }
}


template<typename Set, typename Graph>
void LexBFS::explore(Set &S, Graph &G) {
    initialize(G.size(), S);
    clear();
    
    for(auto pivot : ordering) {
        
#ifdef DBG_LEXBFS
        std::cout << *this << std::endl;
        std::cout << "remove pivot = " << pivot << std::endl;
#endif
        
        auto& pivot_info{vertices_info[pivot]};
        pop_front(pivot_info.interval);
        pivot_info.interval = no_interval;
        
#ifdef DBG_LEXBFS
        std::cout << *this << std::endl;
#endif
        
        for(auto neighbor : G[pivot]) {
            
#ifdef DBG_LEXBFS
            std::cout << " - neighbor " << neighbor ; //<< std::endl;
#endif
            
            auto neighbor_info{vertices_info[neighbor]};
            
//            assert((neighbor_info.interval == no_interval) == (neighbor_info.position < pivot_info.position));
            
//            if(neighbor_info.position < pivot_info.position) {
            if(neighbor_info.interval == no_interval) {
                
#ifdef DBG_LEXBFS
                std::cout << " (already reached)\n" ;
#endif
                
                continue;
            }
            
            auto neighbor_interval{neighbor_info.interval};
            if(intervals_info[neighbor_interval].stamp != pivot) {
                create_interval(neighbor_interval);
                intervals_info[neighbor_interval].stamp = pivot;
            }
    
            move(neighbor);
            
#ifdef DBG_LEXBFS
                std::cout << ": move\n" ;
            std::cout << *this << std::endl;
#endif
        }
    }
}

//}


#endif // _TEMPO_LexBFS_HPP
