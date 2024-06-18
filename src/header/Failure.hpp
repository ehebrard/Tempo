
#ifndef __TEMPO_FAILURE_HPP
#define __TEMPO_FAILURE_HPP

#include <exception>

#include "Constant.hpp"
#include "Explanation.hpp"

namespace tempo {

/**********************************************
* Failure
**********************************************/

template <typename T> class Failure : public std::exception {
public:
  Explanation<T> reason;

  Failure(Explanation<T> r) : reason(r) {}

  virtual const char *what() const throw() { return "Inconsistency (literal)"; }
};

class SearchExhausted : public std::exception {
    
public:
    
    SearchExhausted() = default;
    
  virtual const char *what() const throw() {
    return "Complete search tree exhausted";
  }
};
}

#endif // __FAILURE_HPP
