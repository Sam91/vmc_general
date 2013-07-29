#include "twoflavor.h"

#ifndef U1TWO_H
#define U1TWO_H

//This wave function implements a U(1) algebraic spin liquid for TWO flavors of fermions.
//Here, we suppose an U(1)^2 symmetry, i.e. the flavors are separately conserved (no hybridization in the MF Hamiltonian)
//General hopping matrix for the Fermions can be specified; filling of the two flavors must be identical.

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1;
  double t2; //t' - next nearest hopping

  //hopping matrix between sites per flavor
  double ***t;

  bool ap;  //anti-periodic bc in one direction
};


class u1two : public twoflavor
{
public:
  u1two(int, int);
  ~u1two();

  virtual void print();
  virtual int insert_db();

  void set_hopping(double, double);
  void create_ad();

  parameters* pars;

private:

};

#endif

