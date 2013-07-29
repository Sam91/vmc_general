#include <iostream>

#include "he_two.h"
#include "vmc.h"

//scan over astrow factors that are relevant for a columnar state

int main(int argc, char *argv[])
{
  int req_args = 8;

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
  //wf->set_spiral( (double)atoi(argv[5])/(double)atoi(argv[6]), (double)atoi(argv[7])/(double)atoi(argv[8]) );
  wf->set_spiral( (double)atoi(argv[2])/(double)atoi(argv[3]), (double)atoi(argv[4])/(double)atoi(argv[5]) );

  //wf->print_d();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( 200 ); //number of bins to average over
  
  wf->js[0] = atoi(argv[6 ])/100.; //j1a
  wf->js[1] = atoi(argv[7 ])/100.; //j2a
  wf->js[2] = atoi(argv[8 ])/100.; //j3
  wf->js[3] = atoi(argv[9 ])/100.; //j1b
  wf->js[4] = atoi(argv[10])/100.; //j2b
  wf->js[5] = atoi(argv[11])/100.; //j1c
  wf->js[6] = atoi(argv[12])/100.; //j2c

  //fix isotropic Jastrows
  //wf->js[5] = wf->js[3] = wf->js[0];
  //wf->js[6] = wf->js[4] = wf->js[1];

  for(int n1=0; n1<15; n1++) //sweep over j1_f=j1a, j1_a=j1b=j1c
  {
    wf->js[5] = wf->js[3] = atoi(argv[8])/100.;
    for(int n2=0; n2<15; n2++)
    {
      try
      {
        myvmc->run();
        myvmc->calculate_statistics();
        wf->insert_db();
      } catch (int e) {
        continue;
      }
      wf->js[3] += .05;
      wf->js[5] = wf->js[3];
    }
    wf->js[0] += .05;
  }

  delete myvmc;
  delete wf;

  return 0;
}

