#include "u1hybrid.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]
//Functions/paramaters are defined here that are appropriate for three-sublattices ordering in NS=3 models

#ifndef U1DIRAC_H
#define U1DIRAC_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  bool ap; // define boundary conditions
};


class u1dirac : public u1hybrid
{
public:
  u1dirac(int, int);
  ~u1dirac();

  virtual void print();
  virtual int insert_db();
  virtual void print_avgs();

  void set_hopping();
  void set_hoppingk();
  void set_hopping3(double*, double*);

  void find_max_conf();

  parameters* pars;

private:

};

#endif
