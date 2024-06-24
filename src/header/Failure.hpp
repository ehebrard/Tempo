/************************************************
 * Tempo Failure.hpp
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

#ifndef __TEMPO_FAILURE_HPP
#define __TEMPO_FAILURE_HPP

#include <exception>

#include "Constant.hpp"
#include "Explanation.hpp"

namespace tempo {

//! Failure exception
template <typename T> class Failure : public std::exception {
public:
  Explanation<T> reason;

  Failure(Explanation<T> r) : reason(r) {}

  virtual const char *what() const throw() { return "Inconsistency (literal)"; }
};

//! End-of-search exception
class SearchExhausted : public std::exception {
    
public:
    
    SearchExhausted() = default;
    
  virtual const char *what() const throw() {
    return "Complete search tree exhausted";
  }
};
}

#endif // __FAILURE_HPP
