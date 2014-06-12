#include "wavefunction.h"

#ifndef VMC_H
#define VMC_H

class vmc
{
public:
  ~vmc();
  vmc();

  void set_wf( wavefunction* );
  void thermalize();

  int accepted();

  void initialize(int);

  //acumulates the observables
  void accumulate();

  void run(); //several walks with data collection inbetween

  void calculate_statistics();

private:
  wavefunction *mywf;

  int L, LD, N;
  int nk, runs;

  double **f;
};

#endif

