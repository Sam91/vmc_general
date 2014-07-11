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

  void set_spiral(double, double); //spiral state
  void set_four(double, double); //four sublattice distortion of 120' state

  void set_q0(); //set the q=0 on the kagome

  parameters* pars;

  virtual void print();
//  void find_max_conf();
  virtual int insert_db();
  
//  double jastrow();
//  double jastrow(int, int);
  
protected:

};

#endif
