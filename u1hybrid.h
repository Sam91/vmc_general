#include "wavefunction.h"

//This (abstract) class implements a fermionic unpaired wavefunction for N flavors of fermions. The fermions are hybridized.
//It may be used to test ordering instabillities of the U(1) spin liquid, but the code is more general that that.
//A general hybridization hopping matrix "t" for the Fermions can be specified.
//An ordering patter can be specified with "h" (strength) the on-site projector matrix "d".
//The wave function is projected to a given flavor number; it is appropriate for studying, e.g., SU(N) models.

#ifndef U1HYBRID_H
#define U1HYBRID_H

class u1hybrid : public wavefunction
{
public:
  u1hybrid(int);
  ~u1hybrid();

  virtual void getwf();
  virtual void get_dwf();

#if WFC
  virtual complex<double> swap(int i1, int i2, bool);
  virtual complex<double> crop(int i1, int i2, bool);
#else
  virtual double swap(int i1, int i2, bool);
  virtual double crop(int i1, int i2, bool);
#endif

  void create_ad();

  virtual void backup_data();
  virtual void restore_data();

  void find_max_conf();
  int findmu();
  int findmu(int);
  void set_h(double);
  
  virtual void correct_cff(bool);
  void normalize(double);

  void set_excit( int, int );
private:

#if WFC
  complex<double> **adx, ***dadp, ***dad;
  complex<double> *current_inv, *old_inv;
#else
  double **adx, ***dadp, ***dad;
  double *current_inv, *old_inv;
#endif
  void create_h0();
  void diagonalize();
  void construct_gs();
  void construct_ex();
  void construct_ex1();
  void construct_ex2();
  void construct_ex3();

  int *current_x, *old_x;

  double *N0; //flavor numbers before projection
  
protected:

#if WFC
  complex<double> ****t; //hopping matrix between sites and flavors
  complex<double> **d; //on-site projector to a product state
#else
  double ****t;
  double **d;
#endif

  void set_hopping(double, double, bool);
  void set_hopping3(double*, double*, bool);

  double h;   // H = T - h D (strength of hopping vs projection)
  double *mu; // chemical potentials (d and mu are actually redundant parametrizations to t; but that's ok)

  //mean-field spectrum and eigenvectors
  double *e; //eigenvalues
#if WFC
  complex<double> **v; //eigenvectors
#else
  double **v;
#endif
  double **Na0; //avg flavor number of each eigenvector
  int *occ; //indices of the occupied MF eigenvectors
  bool *isocc;

  int excit;
  int sp; //spin of the excitation
};

#endif

