#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]

//A simple (rotation invariant) real hopping on any lattice

#ifndef U1REAL_H
#define U1REAL_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  //first, second, and third neighbor mean fields
  double *xi; //hoppping amplitude

};


class u1real : public u1hybrid
{
public:
  u1real(int);
  ~u1real();

  virtual void print();
  virtual int insert_db();
  //virtual int insert_file(const char*);
  virtual int insert_file(string);
  virtual void print_avgs();

  void set_hopping();
  void set_hoppingk( double );

  void find_max_conf();

  parameters* pars;

private:

};

#endif
