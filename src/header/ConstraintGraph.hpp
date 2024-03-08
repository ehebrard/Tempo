#ifndef _SCHEDCL_CONSTRAINT_HPP
#define _SCHEDCL_CONSTRAINT_HPP


#include "Graph.hpp"
#include "ResourceConstraint.hpp"

namespace schedcl {



template <typename T> class ConstraintGraph {


public:
  // all constraints
  std::vector<ResourceConstraint<T> *> constraints;
	
	// bipartite graph variables/constraints. successors are triggers, predecessors are scopes
	Graph network;
	
  ConstraintGraph();

  const std::vector<int> &operator[](const int var_id) const;
  // std::vector<int>& scope(const int cons_id);

  // void resize(const size_t n);

  // void initialise_triggers();

  // add a new constraint to the queue, returns the index of that constraint in
  // the queue
  void add(ResourceConstraint<T> *cons);
	int addVar();
	
  void createTrigger(const var v, const int c);
  // void createEventTrigger(const int evt, const int c, const int h);

  // notifies that variable 'id' has changed, activate the corresponding
  // triggers
  void merge(const int x, const int y);

  size_t size() const;
	size_t numVar() const;
	size_t numCon() const;

        bool has(const int u) const;

        //
        // int var(const int i) const;

        std::ostream &display(std::ostream &os) const;
};


template <typename T>
ConstraintGraph<T>::ConstraintGraph() {}

template <typename T>
const std::vector<int> &ConstraintGraph<T>::operator[](const int var_id) const {
  return network[var_id];
}

// template<typename T>
// std::vector<int>& ConstraintGraph<T>::scope(const int cons_id) {
// 	return network[PREDECESSOR][cons_id];
// }

template <typename T> bool ConstraintGraph<T>::has(const int u) const {
  return network.has(u);
}

template<typename T> 
size_t ConstraintGraph<T>::size() const {
	return network.size();
}

template<typename T> 
size_t ConstraintGraph<T>::numCon() const {
	return constraints.size();
}

template<typename T> 
size_t ConstraintGraph<T>::numVar() const {
	return size() - numCon();
}

// add a new variable
template <typename T>
int ConstraintGraph<T>::addVar() {
	network.resize(size()+1);
	return static_cast<int>(size())-1;
}

// add a new constraint 
template <typename T>
void ConstraintGraph<T>::add(ResourceConstraint<T> *cons) {
	constraints.push_back(cons);
	network.resize(numCon()+1);
}

template <typename T>
void ConstraintGraph<T>::createTrigger(const var v, const int c) {
  network.addArc(v, c);
}

template <typename T> void ConstraintGraph<T>::merge(const int x, const int y) {
  network.merge(x, y);
}

template <typename T> std::ostream &ConstraintGraph<T>::display(std::ostream &os) const {

  auto s{numCon()};
  for (auto x{s}; x < s + numVar(); ++x) {
    if (network.has(x) and network.outdegree(x) > 0) {
      std::cout << x << ":";
      for (auto c : network[x]) {
        std::cout << " " << *(constraints[c]);
      }
      std::cout << std::endl;
    }
  }
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const ConstraintGraph<T> &x) {
  return x.display(os);
}
//
// template <class T> std::ostream &operator<<(std::ostream &os, const ResourceConstraint<T> &x) {
//   return x.display(os);
// }

} // namespace SCHEDCL

#endif // _SCHEDCL_SCHEDULING_HPP
