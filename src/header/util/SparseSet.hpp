
#ifndef _TEMPO_SPARSESET_HPP
#define _TEMPO_SPARSESET_HPP

#include <iostream>
#include <vector>

#include "ReversibleObject.hpp"

namespace tempo {


/**********************************************
 * SparseSet
 **********************************************/
/// Sparse set representation

template<typename E=int, typename T=size_t>
class SparseSet {
   
private:
    
    // so that we can remove from the end
    T end_;
    // so that we can remove from the start
    T start_;
    
private:
  /*!@name Parameters*/
  //@{
  /// list of values
  std::vector<E> list_;

  /// values' indices
  std::vector<size_t> index_;
  //@}

public:
  /*!@name Constructors*/
  //@{
  explicit SparseSet(const size_t n = 0);
  explicit SparseSet(const size_t n, BacktrackEnvironment *e);
  //    explicit SparseSet(const T s, const T e, const size_t n=0);
  explicit SparseSet(std::vector<size_t> &common);

  void reserve(const size_t n);

  /*!@name Accessors*/
  //@{
  bool safe_has(const E elt) const;
  bool has(const E elt) const;
	bool isback(const E elt) const;
	bool isfront(const E elt) const;

  size_t capacity() const;
//	size_t count() const;
//	size_t start() const;
  size_t size() const;
    size_t frontsize() const;
    size_t backsize() const;
    size_t start_idx() const;
    size_t end_idx() const;
  bool empty() const;

  void setStart(const T &s);
  void setEnd(const T &s);

  E next(const E elt) const;

  E prev(const E elt) const;

  E operator[](const size_t idx) const;

  E &operator[](const size_t idx);
  //@}

  /*!@name List Manipulation*/
  //@{
  std::vector<E>::iterator begin();
  std::vector<E>::reverse_iterator rbegin();

  std::vector<E>::iterator end();
  std::vector<E>::reverse_iterator rend();

  std::vector<E>::const_iterator begin() const;
  std::vector<E>::const_reverse_iterator rbegin() const;

  std::vector<E>::const_iterator end() const;
  std::vector<E>::const_reverse_iterator rend() const;

  std::vector<E>::iterator fbegin();
  std::vector<E>::reverse_iterator frbegin();

  std::vector<E>::iterator fend();
  std::vector<E>::reverse_iterator frend();

  std::vector<E>::const_iterator fbegin() const;
  std::vector<E>::const_reverse_iterator frbegin() const;

  std::vector<E>::const_iterator fend() const;
  std::vector<E>::const_reverse_iterator frend() const;

  std::vector<E>::iterator bbegin();
  std::vector<E>::reverse_iterator brbegin();

  std::vector<E>::iterator bend();
  std::vector<E>::reverse_iterator brend();

  std::vector<E>::const_iterator bbegin() const;
  std::vector<E>::const_reverse_iterator brbegin() const;

  std::vector<E>::const_iterator bend() const;
  std::vector<E>::const_reverse_iterator brend() const;

  std::vector<E>::iterator get_iterator(const size_t i);
  std::vector<E>::const_iterator get_iterator(const size_t i) const;


  void fill();
	void fill_back();
	void fill_front();

  void clear();

  void resize(const size_t n);

  void pop_back();

  void pop_front();

  void push_back();

  void push_front();

  template <typename Iter>
  void re_index(const Iter beg, const Iter end); // re_index a sublist

  E front() const;

  E back() const;

  E any(const size_t limit=INFTY) const {
    auto m{std::min(limit, size())};
    return list_[start_ + (random() % m)];
  }

  void push_front(const E elt);
	void push_back(const E elt);
  void add(const E elt);
  void safe_add(const E elt);

  void pull_back(const E elt);
  void remove_back(const E elt);
  void front_to_back(const E elt);
  // void safe_remove(const E elt);
  void pull_front(const E elt);
  void remove_front(const E elt);
  void back_to_front(const E elt);

  size_t index(const E elt) const;

  void save_start(size_t &);
  void save_end(size_t &);
  void restore_start(const size_t);
  void restore_end(const size_t);
  //@}
	
  /*!@name Miscellaneous*/
  //@{
  std::ostream &display(std::ostream &os) const;
};


template<typename E, typename T>
SparseSet<E,T>::SparseSet(const size_t n) {
  end_ = 0;
  start_ = 0;
  reserve(n);
}

//
// template<typename E, typename T>
// SparseSet<E,T>::SparseSet(T& s, T& e) : start_(s), end_(e) {
//}

// template<typename E, typename T>
// SparseSet<E,T>::SparseSet(const T s, const T e, const size_t n) : end_(e) ,
// start_(s){
//     reserve(n);
// }

template <typename E, typename T>
SparseSet<E, T>::SparseSet(const size_t n, BacktrackEnvironment *e)
    : end_(0, e), start_(0, e) {
  reserve(n);
}

template<typename E, typename T>
void SparseSet<E,T>::reserve(const size_t n) {
  while (list_.size() < n) {
    index_.push_back(list_.size());
    list_.push_back(list_.size());
  }
}

template<typename E, typename T>
void SparseSet<E,T>::resize(const size_t n) {
  reserve(n);
	fill();
}

//
// void SparseSet<E,T>::save(size_t &stamp1, size_t &stamp2) { stamp1 = end_; stamp2
// = start_; }
// void SparseSet<E,T>::restore(const size_t stamp1, const size_t stamp2) { end_ =
// stamp1; start_ = stamp2; }

template<typename E, typename T>
void SparseSet<E,T>::save_start(size_t &stamp) { stamp = start_; }

template<typename E, typename T>
void SparseSet<E,T>::save_end(size_t &stamp) { stamp = end_; }

template<typename E, typename T>
void SparseSet<E,T>::restore_start(const size_t stamp) { start_ = stamp; }

template<typename E, typename T>
void SparseSet<E,T>::restore_end(const size_t stamp) { end_ = stamp; }

//@}

/*!@name Accessors*/
//@{

template<typename E, typename T>
bool SparseSet<E,T>::safe_has(const E elt) const {
  if (elt >= 0 && (size_t)elt < index_.size())
    return has(elt);
  return false;
}

template <typename E, typename T> bool SparseSet<E,T>::has(const E elt) const {
  return index_[elt] < static_cast<size_t>(end_) and
         index_[elt] >= static_cast<size_t>(start_);
}

template <typename E, typename T> bool SparseSet<E,T>::isfront(const E elt) const {
  return index_[elt] < static_cast<size_t>(start_);
}

template <typename E, typename T> bool SparseSet<E,T>::isback(const E elt) const {
  return index_[elt] >= static_cast<size_t>(end_);
}

template<typename E, typename T>
size_t SparseSet<E,T>::size() const { return end_ - start_; }

template<typename E, typename T>
size_t SparseSet<E,T>::frontsize() const { return start_; }

template<typename E, typename T>
size_t SparseSet<E,T>::backsize() const { return capacity() - end_; }

template<typename E, typename T>
size_t SparseSet<E,T>::start_idx() const { return start_; }

template<typename E, typename T>
size_t SparseSet<E,T>::end_idx() const { return end_; }


//template<typename E, typename T>
//size_t SparseSet<E,T>::last() const { return end_; }
//
//template<typename E, typename T>
//size_t SparseSet<E,T>::start() const { return start_; }

template <typename E, typename T> void SparseSet<E,T>::setStart(const T &s) { start_ = s; }

template <typename E, typename T> void SparseSet<E,T>::setEnd(const T &s) { end_ = s; }

template<typename E, typename T>
size_t SparseSet<E,T>::capacity() const { return index_.size(); }

template<typename E, typename T>
bool SparseSet<E,T>::empty() const { return end_ == start_; }

template<typename E, typename T>
E SparseSet<E,T>::next(const E elt) const {
  size_t idx = index_[elt] + 1;
  return (idx < end_ ? list_[idx] : elt);
}

template<typename E, typename T>
E SparseSet<E,T>::prev(const E elt) const {
  size_t idx = index_[elt];
  return (idx > start_ ? list_[idx - 1] : elt);
}

template<typename E, typename T>
E SparseSet<E,T>::operator[](const size_t idx) const { return list_[idx]; }

template<typename E, typename T>
E &SparseSet<E,T>::operator[](const size_t idx) { return list_[idx]; }
//@}

/*!@name List Manipulation*/
//@{
template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::fbegin() { return list_.begin(); }

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::begin() { return list_.begin() + start_; }

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::bbegin() { return list_.begin() + end_; }

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::frbegin() {
  return list_.rend() - start_;
}

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::rbegin() {
  return list_.rend() - end_;
}

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::brbegin() {
  return list_.rbegin();
}

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::fend() { return list_.begin() + start_; }

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::end() { return list_.begin() + end_; }

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::bend() { return list_.end(); }

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::frend() { return list_.rend(); }

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::rend() {
  return list_.rend() - start_;
}

template<typename E, typename T>
std::vector<E>::reverse_iterator SparseSet<E,T>::brend() {
  return list_.rend() - end_;
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::fbegin() const {
  return list_.begin();
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::begin() const {
  return list_.begin() + start_;
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::bbegin() const {
  return list_.begin() + end_;
}

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::frbegin() const {
  return list_.rend() - start_;
}

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::rbegin() const {
  return list_.rend() - end_;
}

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::brbegin() const {
  return list_.rbegin();
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::fend() const {
  return list_.begin() + start_;
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::end() const {
  return list_.begin() + end_;
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::bend() const { return list_.end(); }

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::frend() const {
  return list_.rend();
}

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::rend() const {
  return list_.rend() - start_;
}

template<typename E, typename T>
std::vector<E>::const_reverse_iterator SparseSet<E,T>::brend() const {
  return list_.rend() - end_;
}

template<typename E, typename T>
std::vector<E>::const_iterator SparseSet<E,T>::get_iterator(const size_t i) const {
	return list_.begin() + i;
}

template<typename E, typename T>
std::vector<E>::iterator SparseSet<E,T>::get_iterator(const size_t i) {
	return list_.begin() + i;
}

// std::vector<E>::iterator SparseSet<E,T>::begin_after() { return end(); }
// std::vector<E>::reverse_iterator SparseSet<E,T>::rbegin_after() {
//   return list_.rend();
// }
//
// std::vector<E>::iterator SparseSet<E,T>::end_after() { return list_.end(); }
// std::vector<E>::reverse_iterator SparseSet<E,T>::rend_after() { return rbegin(); }
//
// std::vector<E>::const_iterator SparseSet<E,T>::begin_after() const {
//   return end();
// }
// std::vector<E>::const_reverse_iterator SparseSet<E,T>::rbegin_after() const {
//   return list_.rend();
// }
//
// std::vector<E>::const_iterator SparseSet<E,T>::end_after() const {
//   return list_.end();
// }
// std::vector<E>::const_reverse_iterator SparseSet<E,T>::rend_after() const {
//   return rend();
// }

template<typename E, typename T>
void SparseSet<E,T>::fill_back() { end_ = list_.size(); }

template<typename E, typename T>
void SparseSet<E,T>::fill_front() { start_ = 0; }

template<typename E, typename T>
void SparseSet<E,T>::fill() { fill_back(); fill_front(); } //size_ = list_.size(); start_ = 0; }

template<typename E, typename T>
void SparseSet<E,T>::clear() { end_ = 0; start_ = 0; }

template<typename E, typename T>
void SparseSet<E,T>::remove_back(const E elt) {
  if (index_[elt] < end_ and index_[elt] >= start_)
    pull_back(elt);
}

template<typename E, typename T>
void SparseSet<E,T>::pull_back(const E elt) {
  auto last = list_[--end_];
  index_[last] = index_[elt];
  list_[index_[elt]] = last;
  list_[end_] = elt;
  index_[elt] = end_;
}

template <typename E, typename T>
void SparseSet<E, T>::front_to_back(const E elt) {
  pull_back(elt);
  ++start_;
}

// void SparseSet<E,T>::safe_remove_front(const E elt) {
//   if (elt >= 0) {
//     if (static_cast<size_t>(elt) >= list_.size()) {
//       reserve(elt + 1);
//     }
//     remove_front(elt);
//   }
// }

template<typename E, typename T>
void SparseSet<E,T>::remove_front(const E elt) {
  if (index_[elt] < end_ and index_[elt] >= start_)
    pull_front(elt);
}

template<typename E, typename T>
void SparseSet<E,T>::pull_front(const E elt) {

  auto first = list_[start_];
  index_[first] = index_[elt];
  list_[index_[elt]] = first;
  list_[start_] = elt;
  index_[elt] = start_;
	++start_;
}

// ...|y..|x..

template <typename E, typename T>
void SparseSet<E, T>::back_to_front(const E elt) {
  pull_front(elt);
  --end_;
}

template <typename E, typename T>
template <typename Iter>
void SparseSet<E, T>::re_index(const Iter beg, const Iter end) {
  for (auto i{beg}; i != end; ++i) {
    index_[*i] = static_cast<T>(i - fbegin());
  }
}

template<typename E, typename T>
void SparseSet<E,T>::pop_back() { --end_; }

template<typename E, typename T>
void SparseSet<E,T>::pop_front() { ++start_; }

template <typename E, typename T> void SparseSet<E, T>::push_back() { ++end_; }

template <typename E, typename T> void SparseSet<E, T>::push_front() {
  --start_;
}

template<typename E, typename T>
E SparseSet<E,T>::front() const { return list_[start_]; }

template<typename E, typename T>
E SparseSet<E,T>::back() const { return list_[end_ - 1]; }

template<typename E, typename T>
void SparseSet<E,T>::safe_add(const E elt) {
  if (elt >= 0) {
    if (static_cast<size_t>(elt) >= list_.size()) {
      reserve(elt + 1);
    }
    add(elt);
  }
}

template<typename E, typename T>
void SparseSet<E,T>::add(const E elt) {
	
	// std::cout << "add " << elt << " to " << *this << std::endl;

        if (index_[elt] >= static_cast<size_t>(end_))
          push_back(elt);
        else if (index_[elt] < static_cast<size_t>(start_))
          push_front(elt);
}

template<typename E, typename T>
void SparseSet<E,T>::push_back(const E elt) {
  auto next = list_[end_];
  index_[next] = index_[elt];
  list_[index_[elt]] = next;
  index_[elt] = end_;
  list_[end_] = elt;
	++end_;
}

template<typename E, typename T>
void SparseSet<E,T>::push_front(const E elt) {
  auto next = list_[--start_];
  index_[next] = index_[elt];
  list_[index_[elt]] = next;
  index_[elt] = start_;
  list_[start_] = elt;
}

template<typename E, typename T>
size_t SparseSet<E,T>::index(const E elt) const { return index_[elt]; }
//@}

template<typename E, typename T>
std::ostream &SparseSet<E,T>::display(std::ostream &os) const {
    for (size_t i{0}; i < capacity(); ++i) {
    if (i == start_)
      os << " |";
    if (i == end_)
      os << " |";
    os << " " << list_[i];
  }
  if (static_cast<size_t>(end_) == capacity())
    os << " |";
  // os << std::endl;

  // os << "(";
  // for (auto it = begin(); it < end(); ++it) {
  //   os << " " << *it;
  // }
  // os << " )";

  return os;
}

template<typename E, typename T>
std::ostream &operator<<(std::ostream &os, const SparseSet<E,T> &x) {
  return x.display(os);
}

template<typename E, typename T>
std::ostream &operator<<(std::ostream &os, const SparseSet<E,T> *x) {
  return (x ? x->display(os) : os);
}

}


#endif // _TEMPO_SPARSESET_HPP
