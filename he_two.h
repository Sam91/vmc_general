#include "huseelser.h"

#ifndef HE_TWO_H
#define HE_TWO_H

//a structure of parameters that define the wave function
struct parameters : public baseparameters
{
  double phi1, phi2; //spiral angles
};

class he_two : public huseelser
{
public:
  he_two(int);
  ~he_two();

  void set_spiral2(double, double); //spiral state (1 D)
  void set_spiral1(double); //spiral state (2 D)
  void set_four(double, double); //four sublattice distortion of 120' state

  void set_q0(); //set the q=0 on the kagome
  void set_sq3(); //set the sq3 state on the kagome
  void set_cbc1(); //set the cuboc-1 state
  void set_cbc2(); //set the cuboc-2 state

  parameters* pars;

  virtual void print();
//  void find_max_conf();
  virtual int insert_db();
  virtual int insert_file( string );
  
//  double jastrow();
//  double jastrow(int, int);
  
protected:

};

#endif
