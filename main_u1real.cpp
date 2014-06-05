#include <iostream>
#include <iomanip>

#include "u1real.h"

int main(int argc, char *argv[])
{
  int req_args = 10;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1real* wf = new u1real( L, 3 );

  wf->pars->apx = atoi(argv[2]); // P/AP boundary conditions
  wf->pars->apy = atoi(argv[3]);

  wf->pars->e2 = atoi(argv[4])==1 ? true : false ; // unit cell doubling
  wf->pars->TR = atoi(argv[5])==1 ? true : false ; // rotation breaking (staggering of hopping in the hexagon)

  wf->pars->xi[0] = ((double)atoi(argv[6]))/100.; // real hopping parameter on one link
  wf->pars->xi[1] = ((double)atoi(argv[7]))/100.;
  wf->pars->xi[2] = ((double)atoi(argv[8]))/100.;

  wf->set_lattice( "kagome" );
  wf->set_hoppingk( ((double)atoi(argv[9]))/100. );
  wf->set_mc_length( 80 );

  //wf->setmu(-.2427);
  if( wf->findmu()>-1 )
  {
    wf->print();

    vmc* myvmc = new vmc();
    myvmc->set_wf( wf );

    myvmc->initialize( atoi(argv[10]) ); //number of bins to average over
    myvmc->run();
    myvmc->calculate_statistics();
    wf->insert_db();

    delete myvmc;
  }

  cout <<"exiting"<<endl;

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
