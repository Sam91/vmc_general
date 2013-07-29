#include <iostream>
#include <iomanip>

#include "lwave2.h"
//#include "amperean.h"

int main(int argc, char *argv[])
{
  int req_args = 6;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);
  
  lwave2* wf = new lwave2( L*L, L );
  //amperean* wf = new amperean( L*L, L );

  //double* dl = new double[4];
  //dl[4] = 22.;

  //some preparatory configuration of the wave function

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  wf->set_lattice( "triangular" );

  wf->pars->apx = atoi(argv[2]); //boundary condition in x-direction
  wf->pars->apy = atoi(argv[3]);
  wf->pars->t1 = (double)atoi(argv[4])/100.; //nearest-neighbor hopping
  wf->pars->dd = (double)atoi(argv[5])/100.; //pairing amplitude
//  wf->pars->dd0 = (double)atoi(argv[5])/100.; //s-pairing for filled states
  wf->pars->phi1 = (double)atoi(argv[6])/100.; //s-pairing for filled states
  wf->pars->phi2 = (double)atoi(argv[7])/100.; //s-pairing for filled states
  if( argc>=9 ) {
    wf->pars->mu = (double)atoi(argv[8])/100.;
    wf->pars->t1b = (double)atoi(argv[9])/100.;
    wf->pars->t1c = (double)atoi(argv[10])/100.;
  } else {
    wf->pars->mu = .8;
    wf->pars->t1b = wf->pars->t1;
    wf->pars->t1c = wf->pars->t1;
  }
//  wf->pars->r  = atoi(argv[6]); // orientation: 0=long axis, 1=short axis
//  wf->pars->lth  = atoi(argv[7])/100.; //controlling the decay of amperean pairing
//  wf->pars->lr1  = atoi(argv[8])/100.;
//  wf->pars->lr2  = atoi(argv[9])/100.;

  //wf->print();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  wf->set_mc_length( 300 );
  myvmc->set_wf( wf );


  wf->findmu();

  //wf->create_dd();
  //wf->create_ad();

  if( argc>=12 ) myvmc->initialize( atoi(argv[11]) ); //number of bins to average over
  else myvmc->initialize( 30 ); //number of bins to average over

  myvmc->run();

  myvmc->calculate_statistics();

  wf->insert_db();

  delete myvmc;
  delete wf;
  return 0;
}

