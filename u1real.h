#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]

//Real hopping on Kagome lattice

#ifndef U1REAL_H
#define U1REAL_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  //first, second, and third neighbor hopping
  double t1, t2, t3;

  //unit-cell doubling
  bool e2;

  //sign of hopping under rotation
  bool tr;
};


class u1real : public u1hybrid
{
public:
  u1real(int, int);
  ~u1real();

  virtual void print();
  virtual int insert_db();
  virtual void print_avgs();

  void set_hopping();
  void set_hoppingk( double );
  void set_hopping3(double*, double*);

  void find_max_conf();

  parameters* pars;

private:

};

#endif
