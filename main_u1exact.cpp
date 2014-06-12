#include <iostream>
#include <iomanip>

#include "u1real.h"
//#include "u1kagome.h"

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

  u1real* wf = new u1real( L, 1 );
  //u1kagome* wf = new u1kagome( L, 1 );

  wf->pars->ap[0] = atoi(argv[2]); // P/AP boundary conditions
  //wf->pars->apy = atoi(argv[3]);

  wf->pars->xi[0] = ((double)atoi(argv[3]))/100.; // real hopping parameter on one link
  wf->pars->xi[1] = ((double)atoi(argv[4]))/100.;
  wf->pars->xi[2] = ((double)atoi(argv[5]))/100.;

  wf->set_lattice( "chain" );
  //wf->set_lattice( "kagome" );

  wf->set_hoppingk( 0. );
  wf->set_mc_length( 80 );

  wf->pars->desc = "U(1) chain (exact)"; 

  //wf->setmu(-.2427);
  if( wf->findmu()>-1 )
  {
    wf->normalize( ((double)atoi(argv[6]))/100. );
    wf->print();

    //vmc* myvmc = new vmc();
    //myvmc->set_wf( wf );

    //myvmc->initialize( atoi(argv[6]) ); //number of bins to average over

    wf->calculate_exact();
    wf->print_f0();

    //myvmc->run();
    //myvmc->calculate_statistics();
    wf->insert_db();

    //delete myvmc;
  }

  cout <<"exiting"<<endl;

  delete wf;
  return 0;
}

