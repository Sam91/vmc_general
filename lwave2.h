#include "paired2k.h"

#ifndef LWAVE2_H
#define LWAVE2_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1, t1b, t1c; //nearest-neighbor hopping
  double t2; //next-neighbor

  double mu; //chemical potential (same for both flavors)

  int* NF; //this wave function is a flavor number eigenstate, and these are the flavor numbers

  double dd; //maximal paring gap
  double lth, lr1, lr2;
  double phi1, phi2;
};

class lwave2 : public paired2k
{
public:
  lwave2(int, int);
  ~lwave2();

  virtual void create_dd(); //sets delta, xi, mu, ap in paired2k
  virtual int insert_db();
  virtual void print();

  void findmu();
  void findmu(int);

  parameters* pars;

};

#endif

