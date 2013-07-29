#include <iostream>
#include <iomanip>

#include "mysql_wrapper.h"

//#include "u1two.h"
#include "u1.h"
#include "vmc.h"


int main(int argc, char* argv[])
{
  int req_args = 2;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

  int L = atoi(argv[1]);

  u1* wf = new u1(L*L, L);

  wf->pars->ap = atoi(argv[2]); //boundary condition in x-direction

  //some preparatory configuration of the wave function

  //wf->set_lattice( "checkerboard" );
  wf->set_lattice( "triangular" );

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  //myvmc.thermalize(); //this is done in run()

  //wf->printalpha();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";
  //wf->getwf();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";

  //cout << "accepted: " << myvmc.accepted() << "\n";

  myvmc->initialize( 10 ); //number of bins to average over

//  for(double t2=0.; t2<=.5; t2+= .5) {
  double t2=0.;
  wf->set_hopping( 1., t2 );
  wf->create_ad();

  myvmc->run();

  myvmc->calculate_statistics();

  wf->insert_db();
//  }

  delete myvmc;
  delete wf;

  return 0;
}

