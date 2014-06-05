#include "u1real.h"

//A unpaired hybridized fermionic wavefunction for NS flavors of fermions [SU(NS) modesl]

//Complex hopping on Kagome lattice

#ifndef U1COMPL_H
#define U1COMPL_H

class u1compl : public u1real
{
public:
  u1compl(int, int);
  ~u1compl();

  void set_hoppingk( double );

private:
};

#endif
