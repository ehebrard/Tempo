#ifndef TEMPO_RESTART
#define TEMPO_RESTART

#include <iostream>


namespace tempo {


 /*! \class RestartPolicy
    \brief  Interface RestartPolicy

    super class for restart-cutoff sequence generators
  */
  class RestartPolicy {
    
  public:

    unsigned int base;
    
    RestartPolicy(const unsigned int b=100);
    virtual ~RestartPolicy() {}

    virtual void reset(unsigned int& limit) = 0;
    virtual void initialize(unsigned int& limit) = 0;
    
  };


  class NoRestart : public RestartPolicy {
    
  public:
    
    NoRestart();
     virtual ~NoRestart() {}
    
    void reset(unsigned int& limit) override ;
    void initialize(unsigned int& limit) override ;
    
  };


  class Geometric : public RestartPolicy {
    
  public:
    
    unsigned int increment;
    double factor;

    Geometric(const unsigned int b=100, const double f=1.333);
    virtual ~Geometric() {}
    
    void reset(unsigned int& limit) override ;
    void initialize(unsigned int& limit) override ;
    
  };

  class Luby : public RestartPolicy {

  private:
    

	int log2_( const unsigned int v ) const;

    unsigned int luby_seq(const int iter) const;
    
    
  public:
    
    unsigned int iteration;

    Luby(const unsigned int b=100);
    virtual ~Luby() {}
    
    void reset(unsigned int& limit) override ;
    void initialize(unsigned int& limit) override ;
    
  };


template<typename S>
class RestartManager {
  
public:

  
    RestartManager(S& s) : caller(s) {
        
//        std::cout << caller.getOptions().restart_policy << std::endl;
        
        if (caller.getOptions().restart_policy == "luby") {
          impl = new Luby(caller.getOptions().restart_base);
        } else if (caller.getOptions().restart_policy == "geom") {
          impl =
              new Geometric(caller.getOptions().restart_base, caller.getOptions().restart_factor);
        } else {
          impl = new NoRestart();
        }
    }
   ~RestartManager() {
       delete impl;
  }

    void reset() {
        impl->reset(restart_limit);
    }
    
    void initialize() {
        impl->initialize(restart_limit);
        restart_limit += caller.num_fails;
    }
    
    bool limit() {
//        
//        std::cout << caller.num_fails << " >? " << restart_limit << std::endl;
        
        return caller.num_fails > restart_limit;
    }
    
    

private:
    unsigned int restart_limit{static_cast<unsigned int>(-1)};
    RestartPolicy *impl;
    S& caller;
    
};

} // namespace
  
#endif
