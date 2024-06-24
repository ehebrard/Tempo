/************************************************
 * Tempo ReversibleObject.hpp
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

#ifndef _TEMPO_REVOBJECT_HPP
#define _TEMPO_REVOBJECT_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

namespace tempo {

class ReversibleObject;

//! Reference all reversible object
/*!
Manage the trail of every reversible object in that environment
 For multithreading, use multiple environments
*/
class BacktrackEnvironment {

public:
  BacktrackEnvironment();

  // the number of saved states
  [[nodiscard]] int level() const;
  // notify the environment that object 'o''s 'undo()' method needs to be called
  // when reverting to the previous saved state
  void save(ReversibleObject *o);
  // save a new state
  void save();

  // restore to state lvl
  void restore(int lvl);
  // restore to the previous state
  void restore();
  // add object 'o' to the list of object that want to be called once on every
  // restore
  void subscribe(ReversibleObject *o);

  void print() const;

  /**
   * Clears everything so that the scheduler can be re-initialized
   */
  void reset();

protected:
  // list of ptrs to reversible objects to call their "undo()"
  std::vector<ReversibleObject *> trail;

  // stamps[l] is the size of the trail at level l, that is, after l calls to
  // "save()"
  std::vector<size_t> stamps;

  // for reversible object that want to be called on every restore, but only
  // once
  std::vector<ReversibleObject *> subscribers;
};

//! Interface for reversible object
class ReversibleObject {
    
protected:
    BacktrackEnvironment *local_env;

public:
    ReversibleObject(BacktrackEnvironment *e) : local_env(e) {};
    virtual ~ReversibleObject() = default;
  ReversibleObject(const ReversibleObject &) = default;
  ReversibleObject(ReversibleObject &&) noexcept = default;
  ReversibleObject &operator=(const ReversibleObject &) = default;
  ReversibleObject &operator=(ReversibleObject &&) noexcept = default;

  // the backtracking method
  virtual void undo() = 0;
  // calls env->save(this)
  void save();
  // called when environment is saved
  virtual void checkpoint();
    
    void setEnv(BacktrackEnvironment *e) {local_env = e;}

  static BacktrackEnvironment *env;
};

//! primitive type
template<typename T>
class Reversible : public ReversibleObject {

public:
  Reversible(BacktrackEnvironment *e=ReversibleObject::env);
  Reversible(T v, BacktrackEnvironment *e=ReversibleObject::env);

  Reversible(const Reversible<T>&) = default;
  Reversible(Reversible<T> &&) noexcept = default;
  Reversible &operator=(Reversible<T> &&) noexcept = default;
  ~Reversible() override = default;


  void undo() override;

  const T &operator()() const;

  Reversible<T> &operator=(const T &rhs);
  Reversible<T> &operator+=(const T &rhs);
  Reversible<T> &operator-=(const T &rhs);
  Reversible<T> &operator++();
  Reversible<T> &operator--();

  bool operator==(const T &rhs) const;
  bool operator!=(const T &rhs) const;
  bool operator<=(const T &rhs) const;
  bool operator<(const T &rhs) const;
  bool operator>=(const T &rhs) const;
  bool operator>(const T &rhs) const;

  Reversible<T> &operator=(const Reversible<T> &x);
  Reversible<T> &operator+=(const Reversible<T> &x);
  Reversible<T> &operator-=(const Reversible<T> &x);

  bool operator==(const Reversible<T> &rhs) const;
  bool operator!=(const Reversible<T> &rhs) const;
  bool operator<=(const Reversible<T> &rhs) const;
  bool operator<(const Reversible<T> &rhs) const;
  bool operator>=(const Reversible<T> &rhs) const;
  bool operator>(const Reversible<T> &rhs) const;

  int getBacktrackLevel(
      const T &val) const; // returns the lowest level to backtrack to so that
                           // the object is stictly more than val
	
	int operator[](int i) const; // returns the i-th value of the reversible object

  operator T() const { return trail.back(); }
  T prev() const {
    assert(trail.size() > 1);
    return trail[trail.size() - 2];
  }

  std::ostream &display(std::ostream &os) const;

private:
  std::vector<T> trail;
  std::vector<int> lvl;
};

//! ???
template<typename T>
class ReversibleReference : public ReversibleObject {

public:
    ReversibleReference(BacktrackEnvironment *e=ReversibleObject::env);
    ReversibleReference(T& ref, BacktrackEnvironment *e=ReversibleObject::env);

    ReversibleReference(const ReversibleReference&) = delete;
    ReversibleReference(ReversibleReference &&) noexcept = default;
    ReversibleReference &operator=(ReversibleReference &&) noexcept = default;
  ~ReversibleReference() override = default;


  void undo() override;

    ReversibleReference<T> &operator=(const T &rhs);
    ReversibleReference<T> &operator+=(const T &rhs);
    ReversibleReference<T> &operator-=(const T &rhs);
    ReversibleReference<T> &operator++();
    ReversibleReference<T> &operator--();

  int getBacktrackLevel(
      const T &val) const; // returns the lowest level to backtrack to so that
                           // the object is stictly more than val
    
    int operator[](int i) const; // returns the i-th value of the reversible object

  T prev() const {
    assert(trail.size() > 1);
    return trail[trail.size() - 2];
  }

  std::ostream &display(std::ostream &os) const;

private:
  std::vector<T> trail;
  std::vector<int> lvl;
    T& reference;
};

//class ReversibleBool : public ReversibleObject {
//
//public:
//    ReversibleBool(BacktrackEnvironment *e=ReversibleObject::env);
//
//    ReversibleBool(const ReversibleBool&) = delete;
//    ReversibleBool(ReversibleBool &&) noexcept = default;
//    ReversibleBool &operator=(ReversibleBool &&) noexcept = default;
//  ~ReversibleBool() override = default;
//
//
//  void undo() override;
//
//  ReversibleBool &operator=(const boolean_state &rhs);
//
//    boolean_state val() const;
//
//  std::ostream &display(std::ostream &os) const;
//
//private:
//    boolean_state value{Unknown};
//};

//! Reversible vector
template <typename T>
class ReversibleVector : public ReversibleObject, public std::vector<T> {

public:
  ReversibleVector(BacktrackEnvironment *e=ReversibleObject::env);
    
    ReversibleVector(std::initializer_list<T> init, BacktrackEnvironment *e=ReversibleObject::env);

  void undo() override;

  void checkpoint() override;

  std::vector<size_t> trail;
};

template <typename T> ReversibleVector<T>::ReversibleVector(BacktrackEnvironment *e) : ReversibleObject(e) {
    local_env->subscribe(this);
}

template <typename T> ReversibleVector<T>::ReversibleVector(std::initializer_list<T> init, BacktrackEnvironment *e) : ReversibleObject(e), std::vector<T>(init) {
    local_env->subscribe(this);
}

template <typename T> void ReversibleVector<T>::undo() {
  std::vector<T>::resize(trail.back());
  trail.pop_back();
}

template <typename T> void ReversibleVector<T>::checkpoint() {
  trail.push_back(std::vector<T>::size());
}



template<typename T>
ReversibleReference<T>::ReversibleReference(T& ref, BacktrackEnvironment *e) : ReversibleObject(e), reference(ref) {
    lvl.push_back(local_env->level());
}

template<typename T>
void ReversibleReference<T>::undo() {
    reference = trail.back();
  trail.pop_back();
  lvl.pop_back();
}

template<typename T>
ReversibleReference<T>& ReversibleReference<T>::operator=(const T& rhs) {
  if (reference != rhs) {
      auto l{local_env->level()};
    if (l > lvl.back()) {
      ReversibleObject::save();
        trail.push_back(reference);
      reference = rhs;
      lvl.push_back(l);
    } else {
      reference = rhs;
    }
  }
  return *this;
}

template<typename T>
ReversibleReference<T>& ReversibleReference<T>::operator+=(const T& rhs) {
  return operator=(trail.back() + rhs);
}

template<typename T>
ReversibleReference<T>& ReversibleReference<T>::operator-=(const T& rhs) {
  return operator=(trail.back() - rhs);
}

template<typename T>
ReversibleReference<T>& ReversibleReference<T>::operator++() {
  return operator=(trail.back() + 1);
}

template<typename T>
ReversibleReference<T>& ReversibleReference<T>::operator--() {
  return operator=(trail.back() - 1);
}


////

template<typename T>
Reversible<T>::Reversible(BacktrackEnvironment *e) : ReversibleObject(e) {
  trail.push_back(0);
    lvl.push_back(local_env->level());
}

template<typename T>
Reversible<T>::Reversible(const T v, BacktrackEnvironment *e) : ReversibleObject(e) {
  trail.push_back(v);
    lvl.push_back(local_env->level());
}

template<typename T>
void Reversible<T>::undo() {
  trail.pop_back();
  lvl.pop_back();
}

template<typename T>
const T& Reversible<T>::operator()() const {
  return trail.back();
}

template<typename T>
Reversible<T>& Reversible<T>::operator=(const T& rhs) {
  if (trail.back() != rhs) {
      auto l{local_env->level()};
    if (l > lvl.back()) {
      ReversibleObject::save();
      trail.push_back(rhs);
      lvl.push_back(l);
    } else {
      trail.back() = rhs;
    }
  }
  return *this;
}

template<typename T>
Reversible<T>& Reversible<T>::operator=(const Reversible<T>& x) {
  *this = x();
  return *this;
}

template<typename T>
Reversible<T>& Reversible<T>::operator+=(const T& rhs) {
  return operator=(trail.back() + rhs);
}

template<typename T>
Reversible<T>& Reversible<T>::operator+=(const Reversible<T>& x) {
  T rhs{trail.back() + x()};
  return operator=(rhs);
}

template<typename T>
Reversible<T>& Reversible<T>::operator-=(const T& rhs) {
  return operator=(trail.back() - rhs);
}

template<typename T>
Reversible<T>& Reversible<T>::operator-=(const Reversible<T>& x) {
  T rhs{trail.back() - x()};
  return operator=(rhs);
}

template<typename T>
Reversible<T>& Reversible<T>::operator++() {
  return operator=(trail.back() + 1);
}

template<typename T>
Reversible<T>& Reversible<T>::operator--() {
  return operator=(trail.back() - 1);
}

template <typename T> bool Reversible<T>::operator==(const T &rhs) const {
  return trail.back() == rhs;
}

template <typename T> bool Reversible<T>::operator!=(const T &rhs) const {
  return trail.back() != rhs;
}

template <typename T> bool Reversible<T>::operator<=(const T &rhs) const {
  return trail.back() <= rhs;
}

template <typename T> bool Reversible<T>::operator<(const T &rhs) const {
  return trail.back() < rhs;
}

template <typename T> bool Reversible<T>::operator>=(const T &rhs) const {
  return trail.back() >= rhs;
}

template <typename T> bool Reversible<T>::operator>(const T &rhs) const {
  return trail.back() > rhs;
}

template <typename T>
bool Reversible<T>::operator==(const Reversible<T> &rhs) const {
  return trail.back() == rhs();
}

template <typename T>
bool Reversible<T>::operator!=(const Reversible<T> &rhs) const {
  return trail.back() != rhs();
}

template <typename T>
bool Reversible<T>::operator<=(const Reversible<T> &rhs) const {
  return trail.back() <= rhs();
}

template <typename T>
bool Reversible<T>::operator<(const Reversible<T> &rhs) const {
  return trail.back() < rhs();
}

template <typename T>
bool Reversible<T>::operator>=(const Reversible<T> &rhs) const {
  return trail.back() >= rhs();
}

template <typename T>
bool Reversible<T>::operator>(const Reversible<T> &rhs) const {
  return trail.back() > rhs();
}

template <typename T> int Reversible<T>::getBacktrackLevel(const T &val) const {
  int i{static_cast<int>(trail.size()) - 1};
  while (i >= 0 and val < trail[i]) {
    --i;
  }
  return lvl[i + 1];
}

template <typename T> int Reversible<T>::operator[](int i) const {
  return trail[(trail.size() +i) % trail.size()];
}

template <typename T> std::ostream &Reversible<T>::display(std::ostream &os) const {
  for (auto x : trail) {
    os << " " << std::setw(5) << x;
  }
  os << std::endl;
  for (auto x : lvl) {
    os << " " << std::setw(5) << x;
  }
  os << std::endl;
  return os;
}

//template <typename T> std::ostream &operator<<(std::ostream &os, const Reversible<T> &x) {
//  return x.display(os);
//}

} // namespace tempo

#endif
