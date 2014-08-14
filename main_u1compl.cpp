#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <unistd.h>
#include "u1kagome.h"

// This is meant to be compiled with WFC=1

int main(int argc, char *argv[])
{
  int req_args = 15;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  cout << "PID: " << getpid() << endl;

  char* name = new char[200]; gethostname(name, 200);
  cout << "hostname: "<< name << endl; delete[] name;

  if(argc-1 != req_args) {
    cout << "Error: incorrect number of arguments\n";
    cout << argv[0] << " L apx apy e2 tR gR smu xi1 xi2 xi3 a1 a2 a3 mc_len nbin\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1kagome* wf = new u1kagome( L );

  wf->pars->ap[0] = atoi(argv[2]); // P/AP boundary conditions
  wf->pars->ap[1] = atoi(argv[3]);

  wf->pars->e2 = atoi(argv[4])==1 ? true : false ; // unit cell doubling
  wf->pars->TR = atoi(argv[5])==1 ? true : false ; // rotation breaking (staggering of hopping in the hexagon)
  wf->pars->gR = atoi(argv[6]); // rotation representation gR

  bool search_mu = atoi(argv[7])==1 ? true : false ; // search for the chemical pot or set it to zero

  wf->pars->xi[0] = ((double)atoi(argv[ 8]))/100.; //all xi fixed
  wf->pars->xi[1] = ((double)atoi(argv[ 9]))/100.;
  wf->pars->xi[2] = ((double)atoi(argv[10]))/100.;

  wf->pars->a[0] = ((double)atoi(argv[11]))/600.; //phase of hopping in units of Pi
  wf->pars->a[1] = ((double)atoi(argv[12]))/600.; 
  wf->pars->a[2] = ((double)atoi(argv[13]))/600.; 

  wf->set_lattice( "kagome" );
  wf->set_mc_length( atoi(argv[14]) );

  string str;
  if( wf->pars->TR )
    str = "c 12";
  else
    str = "c 3";

  if( wf->pars->e2 )
    wf->pars->desc = string("U(1) Dirac; ").append(str);
  else
    wf->pars->desc = string("U(1) FS; ").append(str);

  if(abs(wf->pars->xi[0])<1e-5 ) wf->pars->a[0]=0.;
  if(abs(wf->pars->xi[1])<1e-5 ) wf->pars->a[1]=0.;
  if(abs(wf->pars->xi[2])<1e-5 ) wf->pars->a[2]=0.;

  wf->print();
  wf->set_hoppingk( 0. );

  if( search_mu )
    if( wf->findmu()==-1 ) {
      cout << "Cannot continue..." << endl;
      exit(-1);
    }

  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( atoi(argv[15]) ); //number of bins to average over
  myvmc->run();
  myvmc->calculate_statistics();
  wf->insert_db();
  delete myvmc;

  delete wf;
  return 0;
}

