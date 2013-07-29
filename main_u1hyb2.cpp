#include <iostream>

#include "u1hybtwo.h"
#include "vmc.h"

int main(int argc, char *argv[])
{
  int req_args = 8;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

  if(argc<req_args+1) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);
  
  u1hybtwo* wf = new u1hybtwo( L*L, L );

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  wf->set_lattice( "triangular" );

  wf->pars->ap = atoi(argv[2]); //boundary condition in x-direction

  double* t1 = new double[6];
  double* t2 = new double[6];
  t1[0] = atoi(argv[3])/100.;
  t1[1] = atoi(argv[4])/100.;
  t1[2] = atoi(argv[3])/100.;
  t1[3] = atoi(argv[3])/100.;
  t1[4] = atoi(argv[4])/100.;
  t1[5] = atoi(argv[3])/100.;
  t2[0] = t2[1] = t2[2] = 0.;
  t2[3] = t2[4] = t2[5] = 0.;

  //wf->set_hopping( atoi(argv[3])/100., 0. );
  wf->set_hopping3( t1, t2 );

  wf->set_h( atoi(argv[5])/100. );
  wf->set_spiral( (double)atoi(argv[6])/(double)atoi(argv[7]), (double)atoi(argv[8])/(double)atoi(argv[9]) );

  wf->js[0] = atoi(argv[10])/100.; //j1a
  wf->js[1] = 0.; //j2a
  wf->js[2] = 0.; //j3
  wf->js[3] = atoi(argv[11])/100.; //j1b
  wf->js[4] = 0.; //j2b
  wf->js[5] = atoi(argv[12])/100.; //j1c
  wf->js[6] = 0.; //j2c


  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  wf->findmu();

  wf->create_ad();

  myvmc->initialize( 30 ); //number of bins to average over

  myvmc->run();

  myvmc->calculate_statistics();

  wf->insert_db();

  delete[] t1; delete[] t2;
  delete myvmc;
  delete wf;

  return 0;
}

