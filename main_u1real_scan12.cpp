#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <unistd.h>
#include "u1kagome.h"

//Here, we want to scan the ratios r1/r2 in a reasonable range

int main(int argc, char *argv[])
{
  int req_args = 9;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;
  cout << "PID: " << getpid() << endl;

  if(argc-1<req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1kagome* wf = new u1kagome( L, 3 );

  wf->pars->apx = atoi(argv[2]); // P/AP boundary conditions
  wf->pars->apy = atoi(argv[3]);

  wf->pars->e2 = atoi(argv[4])==1 ? true : false ; // unit cell doubling
  wf->pars->TR = atoi(argv[5])==1 ? true : false ; // rotation breaking (staggering of hopping in the hexagon)

//  wf->pars->xi[0] = ((double)atoi(argv[6]))/100.; // real hopping parameter on one link
//  wf->pars->xi[1] = ((double)atoi(argv[7]))/100.;

  int rtot = atoi(argv[6]); //total xi1 + |xi2|
  int sgn = atoi(argv[7]);  //sign of xi2 to scan
  if( rtot<0 ) {
    cout << "ERROR: We need rmax > 0." << endl;
    exit( -1 );
  }

  wf->pars->xi[2] = ((double)atoi(argv[8]))/100.;

  int r0=1;
  if( sgn == -1 ) {
    cout << "Scanning with xi2 = " << -2 << " ... -" << rtot << "." << endl;
    r0 = 2;
  } else {
    cout << "Scanning with xi2 = " << 0 << " ... " << rtot << "." << endl;
    r0 = 0;
  }

  wf->set_lattice( "kagome" );
  wf->set_mc_length( 80 );

  if( wf->pars->e2 )
    wf->pars->desc = "U(1) Dirac";
  else
    wf->pars->desc = "U(1) FS";

  for( int r=r0; r <= rtot; r += 2 )
  {
    wf->pars->xi[0] = (double)(rtot-r)/100.;
    wf->pars->xi[1] = (double)(sgn*r)/100.;

    wf->print();
    wf->set_hoppingk( 0. );

    if( wf->findmu()>-1 )
    {
      vmc* myvmc = new vmc();
      myvmc->set_wf( wf );

      myvmc->initialize( atoi(argv[9]) ); //number of bins to average over
      myvmc->run();
      myvmc->calculate_statistics();
      wf->insert_db();
      delete myvmc;
    } else
    {
      cout << "skipping..." << endl;
    }
  }

  cout << "exiting." << endl;

  delete wf;
  return 0;

}

