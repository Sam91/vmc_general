#include <iostream>
#include <iomanip>

#include "huseelser.h"
#include "u13.cpp"
#include "vmc.h"

int main()
{
  int L = 12;
  
  huseelser* wf = new huseelser(L);

//some preparatory configuration of the wave function

  //set the 120 state with eta = 0;
  //wf->set120( asin(sqrt(2./3.)) );
  wf->set120( M_PI/2. );

  //Jastrow factors
  wf->pars->j1 = 0.1;
  //wf->pars->j2 = 0.1;

  //wf->print();

//create a vmc object and assign the wf

  vmc myvmc;

  myvmc.set_wf( wf );

  myvmc.thermalize();

  wf->printalpha();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";
  wf->getwf();
  cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";

  //cout << "accepted: " << myvmc.accepted() << "\n";

  myvmc.initialize( 200 ); //number of bins to average over

  myvmc.run();

  myvmc.calculate_statistics();

  return 0;
}
