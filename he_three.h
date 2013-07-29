#include "wavefunction.h"

#ifndef HUSEELSER_H
#define HUSEELSER_H

//a structure of parameters that define the wave function
struct parameters : public baseparameters
{

  //double J1, J2; //nn Jastrow factors
};

  void set120(double); //120 magnetic state
  parameters* pars;


class huseelser : public wavefunction
{
public:
  huseelser(int);
  ~huseelser();

  virtual void getwf();
#if WFC
  virtual complex<double> swap(int i1, int i2);
  virtual complex<double> crop(int i1, int i2);
#else
  virtual double swap(int i1, int i2);
  virtual double crop(int i1, int i2);
#endif

  virtual bool step();
  virtual void find_starting_conf();
  virtual void backup_data();
  virtual void restore_data();

  virtual void print();
  void find_max_conf();
  
protected:
#if WFC
  complex<double>** d; //product state definition
#else
  double** d;
#endif

};

#endif

