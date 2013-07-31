#include <iostream>
#include <iomanip>

#include "u1dirac.h"

int main(int argc, char *argv[])
{
  int req_args = 3;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1dirac* wf = new u1dirac( L, 3 );

  wf->pars->apx = atoi(argv[2]);
  wf->pars->apy = atoi(argv[3]);

  wf->set_lattice( "kagome" );
  wf->set_hoppingk();
  wf->set_mc_length(50);

  wf->findmu();
  wf->print();

  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize(10); //number of bins to average over
  myvmc->run();
  myvmc->calculate_statistics();

  cout <<"exiting"<<endl;

  delete myvmc;
  delete wf;
  return 0;

}
/*
  wf->pars->apx = atoi(argv[2]); //boundary condition in x-direction
  wf->pars->apy = atoi(argv[3]);
  
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

//  wf->insert_db();

  delete myvmc;
  delete wf;
  return 0;
}
*/
