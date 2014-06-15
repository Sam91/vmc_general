#include "twoflavor.h"

#ifndef PAIRED2K_H
#define PAIRED2K_H

class paired2k : public twoflavor
{
public:
  paired2k(int);
  ~paired2k();

  virtual void create_ad();

  virtual void create_dd()=0; //sets delta, xi, mu, ap

private:

protected:

  double *xi;
#if WFC
  complex<double>** delta;
#else
  double** delta;
#endif

  //double mu; //chemical potential

  double *N0; //flavor numbers before projection
};

#endif
