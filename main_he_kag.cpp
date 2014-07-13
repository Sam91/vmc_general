#include <iostream>

#include "he_two.h"
#include "vmc.h"

int main(int argc, char *argv[])
{
  int req_args = 5;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

  if(argc<req_args-1)
  {
    cout << "Error: incorrect number of arguments\n";
    exit(-1);
  }

  int L = atoi(argv[1]);
  
  he_two* wf = new he_two( L );

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  //wf->set_lattice( "triangular" );
  wf->set_lattice( "kagome" );

  //wf->set_q0();
  wf->set_cbc1();
  //wf->set_cbc2();

  //wf->print_d();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( atoi(argv[5]) ); //number of bins to average over
  
  wf->js[0] = atoi(argv[2])/100.; //j1
  wf->js[1] = atoi(argv[3])/100.; //j2
  wf->js[2] = atoi(argv[4])/100.; //j3

  myvmc->run();
  myvmc->calculate_statistics();
  //wf->insert_db();

  //wf->calculate_exact();
  //wf->print_f0();

/*
  for(int n2=0; n2<30; n2++)
  {
      try
      {
        myvmc->run();
        myvmc->calculate_statistics();
        wf->insert_db();
      } catch (int e) {
        continue;
      }
      wf->js[0] += .02;
  }
*/
  
  delete myvmc;
  delete wf;

  return 0;
}

