
#include <iomanip>

#include "util/LexBFS.hpp"


/**********************************************
 * LexBFS
 **********************************************/


const interval_t LexBFS::head = 0;
const interval_t LexBFS::tail = 1;
const interval_t LexBFS::default_interval = 2;
const interval_t LexBFS::no_interval = static_cast<interval_t>(-1);
const vertex_t LexBFS::no_vertex = static_cast<vertex_t>(-1);

LexBFS::LexBFS() {
        next_interval.resize(3);
        prev_interval.resize(3);
        intervals_info.resize(3);
}

//LexBFS::LexBFS(const size_t n) {
//    clear();
//    resize(n);
//}

bool LexBFS::before(vertex_t x, vertex_t y) const {
    return vertices_info[x].position > vertices_info[y].position;
}


void LexBFS::clear() {
    auto n{static_cast<pos_t>(ordering.size())};
//    next_interval.resize(3);
//    prev_interval.resize(3);
//    intervals_info.resize(3);
    next_interval[head] = default_interval;
    prev_interval[tail] = default_interval;
    next_interval[default_interval] = tail;
    prev_interval[default_interval] = head;
    intervals_info[head] = {n,no_vertex};
    intervals_info[tail] = {n,no_vertex};
    intervals_info[default_interval] = {0,no_vertex};
    free_indices.clear();
    for(interval_t i{3}; i<intervals_info.size(); ++i) {
        free_indices.push_back(i);
    }
    
}

void LexBFS::remove(interval_t i) {
    next_interval[prev_interval[i]] = next_interval[i];
    prev_interval[next_interval[i]] = prev_interval[i];
    free_indices.push_back(i);
}

void LexBFS::pop_front(interval_t i) {
    if(++intervals_info[i].first == intervals_info[next_interval[i]].first)
        remove(i);
}

// create an empty interval preceding interval i
void LexBFS::create_interval(interval_t i) {
    interval_t j;
    if(free_indices.empty()) {
        j = static_cast<interval_t>(intervals_info.size());
        intervals_info.emplace_back(intervals_info[i].first, no_vertex);
        next_interval.push_back(i);
        prev_interval.push_back(prev_interval[i]);
    } else {
        j = free_indices.back();
        free_indices.pop_back();
        intervals_info[j] = {intervals_info[i].first, no_vertex};
        next_interval[j] = i;
        prev_interval[j] = prev_interval[i];
    }

    next_interval[prev_interval[i]] = j;
    prev_interval[i] = j;
    
//    intervals_info[i].stamp = pivot;
}

// move a vertex to the previous interval
void LexBFS::move(const vertex_t vertex) {
    auto& vertex_info{vertices_info[vertex]};
    auto start{intervals_info[vertex_info.interval].first};
    auto first_in_interval{ordering[start]};
    
    // swap vertex with the first vertex in its interval
    std::swap(ordering[start], ordering[vertex_info.position]);
    std::swap(vertex_info.position, vertices_info[first_in_interval].position);
    
    // shrink the vertex' original interval
    pop_front(vertex_info.interval);
    
    // move it to the previous intervals
    vertex_info.interval = prev_interval[vertex_info.interval];
}

std::ostream &LexBFS::display(std::ostream &os) const {
    
    
    os << "intervals:";
    for(auto i{next_interval[head]}; i!=tail; i=next_interval[i]) {
        os << " " << i << " (" << intervals_info[i].first << ")";
    }
    os << std::endl;
    
    for(pos_t i{0}; i<ordering.size(); ++i) {
        os << std::setw(3) << i << ": " << std::setw(3) << ordering[i] ;
        if(vertices_info[ordering[i]].interval == no_interval) {
            os << " reached";
        } else {
            auto j{vertices_info[ordering[i]].interval};
            os << " in " << vertices_info[ordering[i]].interval
            << " ([" << intervals_info[j].first
            << ".." << intervals_info[next_interval[j]].first-1
            << "]/" ;
            if(intervals_info[j].stamp == no_vertex)
                os << "none";
            else
                os << intervals_info[j].stamp;
            os << ")";
        }
        os << std::endl;
        assert(vertices_info[ordering[i]].position == i);
        assert(vertices_info[ordering[i]].interval == no_interval or (intervals_info[vertices_info[ordering[i]].interval].first <= i and intervals_info[next_interval[vertices_info[ordering[i]].interval]].first > i));
    }
    
    return os;
}

std::ostream &operator<<(std::ostream &os, const LexBFS &x) {
  return x.display(os);
}

