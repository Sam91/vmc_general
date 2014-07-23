#include <iostream>

#include "he_two.h"
#include "vmc.h"

int main(int argc, char *argv[])
{
  int req_args = 6;

  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

  if(argc<req_args+1)
  {
    cout << "Error: incorrect number of arguments (" << argc << ")\n";
    exit(-1);
  }

  int L = atoi(argv[1]);
  
  he_two* wf = new he_two( L );

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  //wf->set_lattice( "triangular" );
  wf->set_lattice( "kagome" );

  int cbc = atoi(argv[2]);
  if( cbc==1 )
    wf->set_cbc1();
  else if( cbc==2 )
    wf->set_cbc2();
  else if( cbc==3 )
    wf->set_sq3();
  else
    wf->set_q0();

  //wf->print_d();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  myvmc->initialize( atoi(argv[6]) ); //number of bins to average over
  
  wf->js[0] = atoi(argv[3])/100.; //j1
  wf->js[1] = atoi(argv[4])/100.; //j2
  wf->js[2] = atoi(argv[5])/100.; //j3

  //optimal pars for L=20
  //double j1[] = {0,0.02,0.03,0.04,0.06,0.08,0.1,0.12,0.14,0.15,0.16,0.18,0.2,0.21,0.23,0.25,0.26,0.27,0.29,0.31,0.32,0.34,0.35,0.37,0.38,0.4,0.4,0.43,0.44,0.45,0.48};
  //double j2[] = {0.44,0.43,0.42,0.4,0.4,0.38,0.37,0.35,0.34,0.32,0.31,0.29,0.28,0.26,0.25,0.23,0.22,0.2,0.19,0.17,0.16,0.14,0.14,0.12,0.1,0.1,0.08,0.07,0.06,0.04,0.04};
  //double j3[] = {0,0.01,0,0.01,0.02,0.03,0.03,0.04,0.05,0.04,0.04,0.05,0.05,0.05,0.05,0.05,0.06,0.04,0.05,0.04,0.05,0.04,0.04,0.04,0.03,0.04,0.03,0.03,0.02,0.02,0.02};

  //optimal pars for L=22
  //double j1[] = {0.005,0.015,0.03,0.05,0.07,0.08,0.095,0.11,0.13,0.145,0.16,0.18,0.2,0.21,0.22,0.245,0.265,0.275,0.295,0.31,0.325,0.34,0.35,0.37,0.385,0.395,0.41,0.425,0.44,0.45,0.465};
  //double j2[] = {0.46,0.43,0.42,0.4,0.395,0.38,0.37,0.355,0.34,0.32,0.305,0.295,0.28,0.265,0.245,0.235,0.225,0.205,0.19,0.18,0.16,0.155,0.125,0.12,0.11,0.09,0.085,0.075,0.06,0.045,0.03};
  //double j3[] = {0.005,0.005,0.005,0.02,0.025,0.025,0.03,0.035,0.04,0.035,0.045,0.045,0.055,0.055,0.045,0.055,0.055,0.05,0.05,0.045,0.04,0.05,0.035,0.04,0.04,0.035,0.04,0.035,0.03,0.025,0.015};

//  for(int i=0; i<31; i++)
//  {
//    wf->js[0] = j1[i];
//    wf->js[1] = j2[i];
//    wf->js[2] = j3[i];
    myvmc->run();
    myvmc->calculate_statistics();
//    wf->insert_file( argv[3] );
//  }

  wf->insert_db();
  //fappend(wf->js, 3, "test.out");


  //wf->calculate_exact();
  //wf->print_f0();
/*

  for(int n2=0; n2<=70; n2++)
  {
      try
      {
        myvmc->run();
        myvmc->calculate_statistics();
        wf->insert_db();
      } catch (int e) {
        continue;
      }
      if( cbc==0 )
        wf->js[ 0 ] += .02;
      else
        wf->js[ 2 ] += .02;
  }
*/  
  delete myvmc;
  delete wf;

  return 0;
}

