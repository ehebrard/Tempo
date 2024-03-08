
#ifndef _TEMPO_BOUNDSYSTEM_HPP
#define _TEMPO_BOUNDSYSTEM_HPP

#include <iomanip>
#include <iostream>
#include <vector>

#include "Failure.hpp"
#include "Explanation.hpp"
#include "ReversibleObject.hpp"

namespace tempo {


//template<typename T> class Scheduler;

// literals x <= d (interpreted as x >= -d if d is negative)
template <typename T> class BoundConstraint {

public:
    BoundConstraint()
      : l(NoLit), distance(INFTY) {}
    BoundConstraint(const lit l, const T d)
      : l(l), distance(d) {}
    BoundConstraint(const BoundConstraint<T>&) = default;
    BoundConstraint(BoundConstraint<T> &&) noexcept = default;
    BoundConstraint &operator=(const BoundConstraint<T> &) = default;
    BoundConstraint &operator=(BoundConstraint<T> &&) noexcept = default;
    

  lit l;
  T distance; // if distance is negative, interpreted as

    BoundConstraint<T> operator~() const;

  static const BoundConstraint<T> none; //{-1,-1,-1};
  // bool isNull() const { return from<0; }

  bool entails(const BoundConstraint<T>& e) const;
  bool contradicts(const BoundConstraint<T> &e) const;
    
    std::ostream &display(std::ostream &os) const;
};



template <typename T>
class BoundSystem : public ReversibleObject, public Explainer {
  
public:
    
    BoundSystem(
                Scheduler<T>& s,
                const size_t n=2);//, BacktrackEnvironment *e=ReversibleObject::env);
    
    ~BoundSystem() override = default;
    
    void resize(const size_t n);
    
    size_t size() const;
    size_t numLiteral() const;
    
    const BoundConstraint<T>& getConstraint(const size_t i) const;
    
//    T lower(const event x) const;
//    
//    T upper(const event x) const;
    
//    const std::vector<T>& lower() const;
//
//    const std::vector<T>& upper() const;
    
//     std::vector<T>::const_iterator lower() const;
//
//     std::vector<T>::const_iterator upper() const;
    
    const std::vector<T>& lower() const;

    const std::vector<T>& upper() const;
    
    T get(const bool s, const event e) const;
    
    bool falsified(const BoundConstraint<T>& c) const;
    bool satisfied(const BoundConstraint<T>& c) const;
    
    bool set(const bool bt, const event x, const T b, Explanation e=Constant::NoReason);
    
    void undo() override;
    
    lit getLiteral(const lit idx) const;
    lit getIndex(const lit l) const; // get the latest bound-change literal for bound l
    lit getPastIndex(const lit l, const lit s) const; // get the most recent bound-change literal older than s for bound l
    lit getImplicant(const BoundConstraint<T>& c) const; // get the oldest bound-change literal that implies c
    lit getStamp(const lit l) const;
    Explanation getExplanation(const lit) const;
    void xplain(const lit, const hint, std::vector<lit> &) override;// {}
    std::ostream &print_reason(std::ostream &os, const hint) const override;// { return os; }
    int getType() const override;// { return BOUNDEXPL; }
    
    std::ostream &display(std::ostream& os) const;
    std::ostream &displayTrail(std::ostream& os) const;
    
private:
    
    Scheduler<T>& sched;
    
    // bound[LOWER]: distance from each event to the origin (=> -lower_bound)
    // bound[UPPER]: distance from the origin to each event (=> upper_bound)
    std::vector<T> bound[2];
    std::vector<lit> bound_lit[2];
//    std::vector<lit> bound_lit2[2];
//    std::vector<T> bound;
//    std::vector<lit> bound_lit;
    
    // bound literals
    std::vector<BoundConstraint<T>> trail;
    std::vector<lit> stamp;
    std::vector<lit> prev;
    
public:
    std::vector<bool> visited;
    
private:
    std::vector<Explanation> reason;
    
    
    long unsigned int num_prunings{0};
  
};




template <typename T>
bool operator==(const BoundConstraint<T> &d1,
                const BoundConstraint<T> &d2) {
//  return d1.x == d2.x and d1.distance == d2.distance;
    return d1.l == d2.l and d1.distance == d2.distance;
}

template <typename T>
const BoundConstraint<T>
BoundConstraint<T>::none = BoundConstraint<T>(-1, -1);

template <typename T>
BoundConstraint<T> BoundConstraint<T>::operator~() const {
//  return {x, -distance - Gap<T>::epsilon()};
    return {NOT(l), -distance - Gap<T>::epsilon()};
}

template <typename T>
bool BoundConstraint<T>::entails(const BoundConstraint<T>& e) const {
//  return e.x == x and
    return e.l == NOT(l) and
    (distance < 0) == (e.distance < 0) and
    distance <= e.distance;
}

template <typename T>
bool BoundConstraint<T>::contradicts(const BoundConstraint<T> &e) const {
//  return e.x == x and e.distance + distance < 0;
    return e.l == NOT(l) and e.distance + distance < 0;
}

template <typename T>
std::ostream &BoundConstraint<T>::display(std::ostream &os) const {
    os << prettyEvent(EVENT(l));
    if(SIGN(l) == UPPER) {
        os << " <= " << distance;
    } else {
        os << " >= " << -distance;
    }
//    
//    if(distance >= 0) {
//        os << " <= " << distance;
//    } else {
//        os << " >= " << -distance;
//    }
    return os;
}



template <typename T>
BoundSystem<T>::BoundSystem(
                            Scheduler<T>& s,
                            const size_t n)
//, BacktrackEnvironment *e)
:  ReversibleObject(&(s.getEnv())), sched(s)
{
    resize(n); // < 2 ? 2 : n);
}

template <typename T>
void BoundSystem<T>::resize(const size_t n) {
    
//    std::cout << "resize(" << n << ")\n";
    if(n==0)
        return;
    
    auto i{size()};
    
    
//    prev.emplace_back(bound_lit[bt][x]);
//    bound_lit[bt][x] = static_cast<lit>(trail.size());
//    
//    stamp.emplace_back(static_cast<lit>(sched.numEdgeLiteral()-1));
//    trail.emplace_back(LIT(x,bt), b); //bound[bt][x]);
//    reason.emplace_back(e);
//    bound[bt][x] = b;
    
    
    
    bound[LOWER].resize(n, 0);
    bound[UPPER].resize(n, INFTY/2);
    
    bound_lit[LOWER].resize(n, -1);
    bound_lit[UPPER].resize(n, -1);
    bound[UPPER][ORIGIN] = 0;
    
//    visited.resize(2*n, false);
    
    while(i < n) {
        
//        std::cout << i << "/" << bound_lit[LOWER].size() << std::endl;
        
        visited.emplace_back(false);
        prev.emplace_back(NoLit);
        bound_lit[LOWER][i] = static_cast<lit>(trail.size());
        trail.emplace_back(LIT(i,LOWER), 0);
        reason.emplace_back(Constant::NoReason);
        stamp.emplace_back(NoLit);
        
        visited.emplace_back(false);
        prev.emplace_back(NoLit);
        bound_lit[UPPER][i] = static_cast<lit>(trail.size());
        trail.emplace_back(LIT(i,UPPER), INFTY/2);
        reason.emplace_back(Constant::NoReason);
        stamp.emplace_back(NoLit);
        
        ++i;
    }
}


template<typename T>
bool BoundSystem<T>::falsified(const BoundConstraint<T>& c) const {
    return bound[NOT(SIGN(c.l))][EVENT(c.l)] + c.distance < 0;
}

template<typename T>
bool BoundSystem<T>::satisfied(const BoundConstraint<T>& c) const {
    return bound[SIGN(c.l)][EVENT(c.l)] <= c.distance;
}

template <typename T>
const BoundConstraint<T>& BoundSystem<T>::getConstraint(const size_t i) const {
    
    assert(i >= 0 and i < trail.size());
    
    return trail[i];
}

template <typename T>
bool BoundSystem<T>::set(const bool bt, const event x, const T b, Explanation e) {

    
    if(bound[bt][x] > b) {
        
        ++num_prunings;
        
#ifdef DBG_TRACE
        if (DBG_TRACE & PROPAGATION) {
            std::cout << "pruning: ";
            if(bt)
                std::cout << prettyEvent(x) << " <= " << b ;
            else
                std::cout << prettyEvent(x) << " >= " << -b ;
            if (DBG_TRACE & LEARNING) {
                std::cout << " b/c " << e;
            }
            std::cout << std::endl;
        }
#endif
        
        ReversibleObject::save();
        
        prev.emplace_back(bound_lit[bt][x]);
        bound_lit[bt][x] = static_cast<lit>(numLiteral());
        
        visited.push_back(false);
        
//        assert(sched.getIndex(e.the_hint) == (sched.numEdgeLiteral()-1));

//
//        if(e.the_hint >= 0 and sched.getIndex(VAR(e.the_hint)) != (sched.numEdgeLiteral()-1)) {
//            std::cout << "there " << e.the_hint << std::endl;
//            std::cout << "here!!! (" << sched.getIndex(VAR(e.the_hint)) << "/" << (sched.numEdgeLiteral()-1) << ")\n";
//            exit(1);
//        }
        
        stamp.emplace_back(static_cast<lit>(sched.numEdgeLiteral()-1));
        trail.emplace_back(LIT(x,bt), b); //bound[bt][x]);
        reason.emplace_back(e);
        bound[bt][x] = b;
        
        assert(visited.size() == stamp.size());
        assert(trail.size() == stamp.size());
        assert(reason.size() == stamp.size());
        
        if(b + bound[NOT(bt)][x] < 0) {
            
#ifdef DBG_TRACE
            if (DBG_TRACE & PROPAGATION) {
                std::cout << "FAIL on bound!\n";
            }
#endif
            
            throw Failure({this, LIT(x, bt)});
            
        }
 
        return true;
    }
    
    return false;
}

template <typename T> void BoundSystem<T>::undo() {
    
    assert(trail.size() == reason.size());
    
//    BoundConstraint<T> b{trail.back()};
    BoundConstraint<T> b{trail[prev.back()]};
    auto bt{SIGN(b.l)};
    auto x{EVENT(b.l)};
    bound[bt][x] = b.distance;
    bound_lit[bt][x] = prev.back();
    
    prev.pop_back();
    trail.pop_back();
    stamp.pop_back();
    reason.pop_back();
    visited.pop_back();
}

template <typename T>
size_t BoundSystem<T>::size() const {
    return bound[LOWER].size();
//    return bound.size()/2;
}

template <typename T>
size_t BoundSystem<T>::numLiteral() const {
    return trail.size();
}

//template <typename T>
// std::vector<T>::const_iterator BoundSystem<T>::lower() const {
//    return bound.begin();
//}
//
//template <typename T>
// std::vector<T>::const_iterator BoundSystem<T>::upper() const {
//    return bound.begin()+size();
//}

template <typename T>
 const std::vector<T>& BoundSystem<T>::lower() const {
     return bound[LOWER];
}

template <typename T>
 const std::vector<T>& BoundSystem<T>::upper() const {
     return bound[UPPER];
}

template <typename T>
T BoundSystem<T>::get(const bool s, const event e) const {
    return bound[s][e];
}

template <typename T>
lit BoundSystem<T>::getLiteral(const lit idx) const {
    return getConstraint(idx).l;
}

template <typename T>
lit BoundSystem<T>::getIndex(const lit l) const {
    return bound_lit[SIGN(l)][EVENT(l)];
}

template <typename T>
lit BoundSystem<T>::getPastIndex(const lit l, const lit s) const {
    
//    std::cout << "\nfind lit for " << (SIGN(l)==LOWER ? "l" : "u") << "b(" << prettyEvent(EVENT(l)) << ") older than " << s << std::endl;
    
    auto p{bound_lit[SIGN(l)][EVENT(l)]};
    
//    std::cout << " - " << trail[p] << "?\n";
    while(p >= s) {
        p = prev[p];
//        std::cout << " - " << trail[p] << "?\n";
    }
    return p;
}

template <typename T>
lit BoundSystem<T>::getImplicant(const BoundConstraint<T>& c) const {
//    std::cout << "\nfind lit implying " << c << std::endl;
    auto p{bound_lit[SIGN(c.l)][EVENT(c.l)]};
    
    auto r{p};
//    std::cout << " - " << trail[p] << "?\n";
    while(p != NoLit and trail[p].distance <= c.distance) {
//        std::cout << " - " << trail[p] << "?\n";
        r = p;
        p = prev[p];
    }

    return r;
}

template <typename T>
lit BoundSystem<T>::getStamp(const lit l) const {
    
    assert(stamp.size() == trail.size());
    
    assert(l < static_cast<lit>(stamp.size()));
    
    return stamp[l];
}

template <typename T>
Explanation BoundSystem<T>::getExplanation(const lit l) const {
    return reason[l];
}

template<typename T>
void BoundSystem<T>::xplain(const lit l, const hint h, std::vector<lit> &Cl) {
    if(l == NoLit) {
        auto bt{SIGN(h)};
        auto x{EVENT(h)};
        Cl.push_back(BOUND(bound_lit[bt][x]));
        Cl.push_back(BOUND(bound_lit[NOT(bt)][x]));
    } else {
        assert(false);
        //std::cout << "bound revision from edge [" << sched.getEdge(h) << "]";
    }
}

template<typename T>
std::ostream &BoundSystem<T>::print_reason(std::ostream &os, const hint) const {
    os << "bounds";
    return os;
}

template<typename T>
int BoundSystem<T>::getType() const { return BOUNDEXPL; }

template <typename T>
std::ostream &BoundSystem<T>::display(std::ostream& os) const {
    
    for(size_t x{2}; x<size(); ++x) {
        os << std::setw(3) << std::left << prettyEvent(x) << ": [" << (lower()[x] ? -lower()[x] : 0) << ".." << upper()[x] << "]\n";
    }
    os << "makespan: [" << -lower()[1] << ".." << upper()[1] << "]\n";
    
    return os;
}

template <typename T>
std::ostream &BoundSystem<T>::displayTrail(std::ostream& os) const {
    
    for(size_t i{0}; i<trail.size(); ++i)
    {
        auto c{trail[i]};
        os << std::setw(5) << i << std::setw(5) << stamp[i] << std::endl;
    }
    return os;
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const BoundSystem<T> &x) {
  return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const BoundConstraint<T> &x) {
  return x.display(os);
}

} // namespace tempo

#endif
