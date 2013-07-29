#include "wavefunction.h"

#ifndef U1_H
#define U1_H

//This wave function implements a U(1) algebraic spin liquid for NS flavors of fermions.
//Here, we suppose an U(1)^N symmetry, i.e. all flavors are separately conserved (no hybridization in the MF Hamiltonian)
//General hopping matrix for the Fermions can be specified

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1;
  double t2; //t' - next nearest hopping

  //hopping matrix between sites per flavor
  double ***t;

  int* NF; //this wave function is a flavor number eigenstate, and these are the flavor numbers

  bool ap;  //anti-periodic bc in one direction
};


class u1 : public wavefunction
{
public:
  u1(int, int);
  ~u1();

  virtual void getwf();
#if WFC
  virtual complex<double> swap(int i1, int i2, bool);
  virtual complex<double> crop(int i1, int i2, bool);
#else
  virtual double swap(int i1, int i2, bool);
  virtual double crop(int i1, int i2, bool);
#endif

  virtual void find_starting_conf();
  virtual void backup_data();
  void backup_data(int, int);
  virtual void restore_data();

  virtual void print();
  virtual int insert_db();

  void set_hopping(double, double);
  void create_ad();
  void find_max_conf();

  parameters* pars;

private:

  double** adx;
  double** current_inv;
  double** old_inv;

  int **current_x;
  int **old_x;

  int bkp_f1, bkp_f2;

#if WFC
  complex<double> *wfn, *wfn_old;
#else
  double *wfn, *wfn_old;
#endif
  double *cff; //wave function normalization
  void normalize(double);
  void normalize(double, int);
  void correct_cff( bool );
};

#endif

