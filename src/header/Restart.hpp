#ifndef TEMPO_RESTART
#define TEMPO_RESTART

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

} // namespace
  
#endif
