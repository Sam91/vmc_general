#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]
//Functions/paramaters are defined here that are appropriate for three-sublattices ordering in NS=3 models

#ifndef U1HYBTWO_H
#define U1HYBTWO_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1, t1b, t1c;
  double t2, t2b, t2c; //t' - next nearest hopping


  double phi1; //may be used to define an ordering pattern (d)
  double phi2; //other ordering pattern angles

  bool ap;
};


class u1hybtwo : public u1hybrid
{
public:
  u1hybtwo(int);
  ~u1hybtwo();

  virtual void print();
  virtual int insert_db();

  void set_hopping(double, double);
  void set_hopping3(double*, double*);
  void set_spiral(double, double);

  void find_max_conf();

  parameters* pars;

private:

};

#endif
