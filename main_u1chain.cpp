#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "u1real.h"

int main(int argc, char *argv[])
{
  int req_args = 8;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  cout << "PID: " << getpid() << endl;

  char* name = new char[200]; gethostname(name, 200);
  cout << "hostname: "<< name << endl; delete[] name;

  if(argc-1 != req_args) {
    cout << "Error: incorrect number of arguments\n";
    cout << argv[0] << " L apx smu xi1 xi2 xi3 mc_len nbin\n";
    exit(-1);
  }

  int L = atoi(argv[1]);

  u1real* wf = new u1real( L );

  wf->pars->ap[0] = atoi(argv[2]); // P/AP boundary conditions

  bool search_mu = false; //atoi(argv[3])==1 ? true : false ; // search for the chemical pot or set it to zero

  wf->pars->xi[0] = ((double)atoi(argv[4]))/100.; //all xi fixed
  wf->pars->xi[1] = ((double)atoi(argv[5]))/100.;
  wf->pars->xi[2] = ((double)atoi(argv[6]))/100.;

  wf->set_lattice( "chain" );

  wf->set_hoppingk( 0. );
  wf->set_mc_length( atoi(argv[7]) );

  wf->pars->desc = "U(1) chain"; 
  wf->set_excit( atoi(argv[3]), 0 );

  if( search_mu ) {
    if( wf->findmu()==-1 ) {
      cout << "Cannot continue..." << endl;
      exit(-1);
    }
  } else
    wf->create_ad();

  wf->print();

  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( atoi(argv[8] ) ); //number of bins to average over
  myvmc->run();
  myvmc->calculate_statistics();

  //char* fn = new char[100];
  string fn = "/users/invites/sbieri/chain_tst.out";

  wf->insert_db();
  //wf->insert_file( fn );

  delete myvmc;

  delete wf;
  return 0;
}

