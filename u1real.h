#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]

//Real hopping on Kagome lattice

#ifndef U1REAL_H
#define U1REAL_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  //first, second, and third neighbor mean fields
  double *xi; //hoppping amplitude
  double *dd; //pairing amplitude

  double *a; //hopping phase
  double *b; //pairing phase

  double *ll; //on-site term (chemical potential)
  
  bool e2; //unit-cell doubling

  bool TR; //sign of MF under rotation
  int gR;  //SU(2) representation of rotation
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
  bool rr; //some local parameter for U(1) states that tells me if the real part of the hopping is staggered
};

#endif
