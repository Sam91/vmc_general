#include <iostream>

#include "he_two.h"
#include "vmc.h"

int main(int argc, char *argv[])
{
  int req_args = 13;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

  if(argc<req_args-1)
  {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);
  
  he_two* wf = new he_two( L*L, L );

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  wf->set_lattice( "triangular" );
    
  //wf->set_four( (double)atoi(argv[5])/(double)atoi(argv[6]), (double)atoi(argv[7])/(double)atoi(argv[8]) );
  wf->set_spiral( (double)atoi(argv[2])/(double)atoi(argv[3]), (double)atoi(argv[4])/(double)atoi(argv[5]) );

  //wf->print_d();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( 200 ); //number of bins to average over
  
  wf->js[0] = atoi(argv[6 ])/100.;
  wf->js[1] = atoi(argv[7 ])/100.;
  wf->js[2] = atoi(argv[8 ])/100.;
  wf->js[3] = atoi(argv[9 ])/100.;
  wf->js[4] = atoi(argv[10])/100.;
  wf->js[5] = atoi(argv[11])/100.;
  wf->js[6] = atoi(argv[12])/100.;

  wf->jl = atoi(argv[13])/100.;

  myvmc->run();
  myvmc->calculate_statistics();
  wf->insert_db();

  delete myvmc;
  delete wf;

  return 0;
}

