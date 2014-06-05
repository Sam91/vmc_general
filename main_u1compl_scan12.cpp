#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <unistd.h>
#include "u1compl.h"

//Here, we want to scant the ratios r1/r2 in a reasonable range

int main(int argc, char *argv[])
{
  int req_args = 11;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;
  cout << "PID: " << getpid() << endl;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1compl* wf = new u1compl( L, 3 );

  wf->pars->apx = atoi(argv[2]); // P/AP boundary conditions
  wf->pars->apy = atoi(argv[3]);

  wf->pars->e2 = atoi(argv[4])==1 ? true : false ; // unit cell doubling
  wf->pars->TR = atoi(argv[5])==1 ? true : false ; // rotation breaking (staggering of hopping in the hexagon)

//  wf->pars->xi[0] = ((double)atoi(argv[6]))/100.; // real hopping parameter on one link
//  wf->pars->xi[1] = ((double)atoi(argv[7]))/100.;
  wf->pars->xi[2] = ((double)atoi(argv[9]))/100.;

  int r0 = atoi(argv[6]);
  int r1 = atoi(argv[7]);
  int rtot = atoi(argv[8]);
  if( r1<r0 ) {
    cout << "ERROR: We need r1 > r0." << endl;
    exit( -1 );
  }

  cout << "Scanning with r = " << r0 << " ... " << r1 << "; rtot=" << rtot << endl;

  wf->set_lattice( "kagome" );
  wf->set_mc_length( 80 );

  if( wf->pars->e2 )
    wf->pars->desc = "U(1) Dirac - c";
  else
    wf->pars->desc = "U(1) FS - c";

  for( int r=r0; r<=r1; r+=2 )
  {
    wf->pars->xi[0] = (double)(rtot-abs(r))/100.;
    wf->pars->xi[1] = (double)r/100.;

    wf->print();
    wf->set_hoppingk( ((double)atoi(argv[10]))/100. );

    if( wf->findmu()>-1 )
    {
      vmc* myvmc = new vmc();
      myvmc->set_wf( wf );

      myvmc->initialize( atoi(argv[11]) ); //number of bins to average over
      myvmc->run();
      myvmc->calculate_statistics();
      wf->insert_db();
      delete myvmc;
    } else
    {
      cout << "skipping..." << endl;
    }
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
