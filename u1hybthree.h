#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]
//Functions/paramaters are defined here that are appropriate for three-sublattices ordering in NS=3 models

#ifndef U1HYBTHREE_H
#define U1HYBTHREE_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1;
  double t2; //t' - next nearest hopping

  double r2; // additional next neighbor hopping on the right link

  double phi1; //may be used to define an ordering pattern (d)
  double eta, phi, theta, psi; //other ordering pattern angles

  double td1, td2; //diagonal (nn) hoppings on the square lattice

  bool ap;
  
//These parameters are now defined in u1hybrid
//  double *mu; //chemical potential of the flavors
//  int* NF; //the wave function is projected to this flavor number
//  double h; //strength of hopping vs on site projector
};


class u1hybthree : public u1hybrid
{
public:
  u1hybthree(int, int);
  ~u1hybthree();

  virtual void print();
  virtual int insert_db();

  void set_hopping(double, double);
  void set_umb(double);
  void set_three(double);
  void set_three(double,double,double); //real three-sublattice state
  void set_three(double, double, double,double,double); //complex three-sublattice state
  void set_four(double);
  void set_xxyz(int);

  void find_max_conf();

  void set_td(double, double);
  void set_r2(double);

  parameters* pars;

private:

};

#endif
