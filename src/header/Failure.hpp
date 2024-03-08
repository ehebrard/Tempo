
#ifndef __TEMPO_FAILURE_HPP
#define __TEMPO_FAILURE_HPP

#include <exception>

#include "Explanation.hpp"

namespace tempo {

/**********************************************
* Failure
**********************************************/

class Failure: public std::exception
{
public:
    
    Explanation reason;
    
    Failure(Explanation r=Constant::NoReason) : reason(r) {}


  virtual const char* what() const throw()
  {
    return "Inconsistency (literal)";
  }
};

//class NegativeCycle: public std::exception
//{
//public:
//    NegativeCycle() = default;
//
//  virtual const char* what() const throw()
//  {
//    return "Inconsistency (negative cycle)";
//  }
//};

class SearchExhausted : public std::exception {
    
public:
    
    SearchExhausted() = default;
    
  virtual const char *what() const throw() {
    return "Complete search tree exhausted";
  }
};
}

#endif // __FAILURE_HPP
