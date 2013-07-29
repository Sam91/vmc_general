#include "wavefunction.h"

#ifndef HUSEELSER_H
#define HUSEELSER_H

class huseelser : public wavefunction
{
public:
  huseelser(int, int);
  ~huseelser();

  virtual void getwf();
#if WFC
  virtual complex<double> swap(int i1, int i2, bool);
  virtual complex<double> crop(int i1, int i2, bool);
#else
  virtual double swap(int i1, int i2, bool);
  virtual double crop(int i1, int i2, bool);
#endif

  virtual void backup_data();
  virtual void restore_data();

  virtual void correct_cff(bool);
  void normalize(double);

  void find_max_conf();
  void print_d();
  
protected:
#if WFC
  complex<double>** d; //product state definition
#else
  double** d;
#endif

};

#endif
