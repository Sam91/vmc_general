#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <unistd.h>
#include "u1kagome.h"

//Here, we want to scan the ratios r1/r3 in a reasonable range
// This is meant to be compiled with WFC=1

int main(int argc, char *argv[])
{
#if WFC==0
  cout << "ERROR: Need to compuile with WFC=1\n";
  exit(-1);
#endif

  int req_args = 14;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;
  cout << "PID: " << getpid() << endl;

  if(argc-1 != req_args) {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1kagome* wf = new u1kagome( L, 3 );

  wf->pars->apx = atoi(argv[2]); // P/AP boundary conditions
  wf->pars->apy = atoi(argv[3]);

  wf->pars->e2 = atoi(argv[4])==1 ? true : false ; // unit cell doubling
  wf->pars->TR = atoi(argv[5])==1 ? true : false ; // rotation breaking (staggering of hopping in the hexagon)

  wf->pars->gR = atoi(argv[6]); // rotation representation gR

//  wf->pars->xi[0] = ((double)atoi(argv[6]))/100.; // real hopping parameter on one link
//  wf->pars->xi[1] = ((double)atoi(argv[7]))/100.;

  bool search_mu = atoi(argv[7])==1 ? true : false ; // search for the chemical pot or set it to zero

  int rtot = atoi(argv[8]); //total xi1 + |xi2|
  int sgn = atoi(argv[9]);  //sign of xi2 to scan
  if( rtot<0 ) {
    cout << "ERROR: We need rmax > 0." << endl;
    exit( -1 );
  }

  wf->pars->xi[2] = ((double)atoi(argv[10]))/100.; //xi2 is fixed

  wf->pars->a[0] = ((double)atoi(argv[11]))/600.; //phase of hopping in units of Pi
  wf->pars->a[1] = ((double)atoi(argv[12]))/600.; 
  wf->pars->a[2] = ((double)atoi(argv[13]))/600.; 

  int r0=1;
  if( sgn == -1 ) {
    cout << "Scanning with xi3 = " << -2 << " ... -" << rtot << "." << endl;
    r0 = 2;
  } else { //start x2 at 2 for now
    cout << "Scanning with xi3 = " << 2 << " ... " << rtot << "." << endl;
    r0 = 2;
  }

  wf->set_lattice( "kagome" );
  wf->set_mc_length( 80 );

  string str;
  if( wf->pars->TR )
    str = "c 12";
  else
    str = "c 3";

  if( wf->pars->e2 )
    wf->pars->desc = string("U(1) Dirac: ").append(str);
  else
    wf->pars->desc = string("U(1) FS; ").append(str);


  for( int r=r0; r <= rtot; r += 2 )
  {
    wf->pars->xi[0] = (double)(rtot-r)/100.;
    wf->pars->xi[2] = (double)(sgn*r)/100.;

    if(abs(wf->pars->xi[0])<1e-5 ) wf->pars->a[0]=0.;
    if(abs(wf->pars->xi[1])<1e-5 ) wf->pars->a[1]=0.;
    if(abs(wf->pars->xi[2])<1e-5 ) wf->pars->a[2]=0.;

    wf->print();
    wf->set_hoppingk( 0. );

    if( search_mu )
    {
      if( wf->findmu()==-1 ) {
        cout << "skipping..." << endl;
        continue;
      }
    } else wf->create_ad();

    vmc* myvmc = new vmc();
    myvmc->set_wf( wf );

    myvmc->initialize( atoi(argv[14]) ); //number of bins to average over
    myvmc->run();
    myvmc->calculate_statistics();
    wf->insert_db();
    delete myvmc;
  }

  cout << "exiting." << endl;

  delete wf;
  return 0;

}

