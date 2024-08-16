

#include "Global.hpp"
#include "Restart.hpp"


using namespace tempo;


    RestartPolicy::RestartPolicy(const unsigned int b) : base(b) {}


    NoRestart::NoRestart() : RestartPolicy(std::numeric_limits<unsigned int>::max()) {}
    
    void NoRestart::reset(unsigned int& limit) {
      limit = base;
    }

    void NoRestart::initialize(unsigned int& limit, const unsigned int) {
      limit = base;
    }
  

    Geometric::Geometric(const unsigned int b, const double f) : RestartPolicy(b), factor(f) {}
    
    void Geometric::reset(unsigned int& limit) {
      limit += increment;
      increment = (unsigned int)((double)increment * factor);
    }

    void Geometric::initialize(unsigned int& limit, const unsigned int val) {
      limit = 0;
      increment = base;
      reset(limit);
        limit += val;
    }


    int Luby::log2_( const unsigned int v ) const {
	  if( !v ) return -std::numeric_limits<int>::max();
	  int exponent = -1;
	  while( (v >> (++exponent)) > 1 );
	  return exponent;
	}

    
    unsigned int Luby::luby_seq(const int iter) const {
      unsigned int thelog = log2_(iter);
      if( iter == (1 << (thelog + 1))-1 )
		return (1 << thelog);
      return luby_seq(iter - (1 << thelog) + 1);
    }


    Luby::Luby(const unsigned int b) : RestartPolicy(b) { iteration=0; }
    
    void Luby::reset(unsigned int& limit) {
      limit += (base * luby_seq(++iteration));
    }

    void Luby::initialize(unsigned int& limit, const unsigned int val) {
      iteration = 0;
      reset(limit);
        limit += val;
    }
 
