#include "wavefunction.h"

#ifndef TWOFLAVOR_H
#define TWOFLAVOR_H

/*
 * This is a base class for implementing wave functions for two flavors with N0=N1.
 *
 * The wave function must be fully paired or unpaired, such that it can be written as a single determinant.
 */


class twoflavor : public wavefunction
{
public:
  twoflavor(int, int);
  ~twoflavor();

  virtual void getwf();
#if WFC
  virtual complex<double> swap(int i1, int i2, bool);
  virtual complex<double> crop(int i1, int i2, bool);
#else
  virtual double swap(int i1, int i2, bool);
  virtual double crop(int i1, int i2, bool);
#endif

  virtual void find_starting_conf();

  void backup_data();
  virtual void restore_data();

  virtual void create_ad() =0;
  virtual void correct_cff(bool);

protected:
#if WFC
  complex<double>* adx;
  complex<double> *current_inv, *old_inv;
  complex<double>* sample_row;
#else
  double* adx; //adx does not have a flavor index here
  double *current_inv, *old_inv;
  double* sample_row;
#endif

  int **current_x, **old_x;

  void normalize(double); // normalize the adx matrix  
  double cff; //wf normalization
  int N2; // N/2, particle number
};

#endif

