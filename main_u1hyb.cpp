#include <iostream>

#include "u1hybtwo.h"
#include "vmc.h"

int main(int argc, char *argv[])
{
  for(int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n";

//if(argc<11) {
//    cout << "Error: incorrect number of arguments\n";
//    exit(-1);
//  }

  int L = atoi(argv[1]);
  
  u1hybthree* wf = new u1hybthree( L*L, L );

  //some preparatory configuration of the wave function

  //wf->set_lattice( "square" );
  //wf->set_lattice( "checkerboard" );
  wf->set_lattice( "triangular" );

  wf->set_h( atoi(argv[5])/100. );
  //wf->set_umb( atoi(argv[5])/100. );
  wf->set_three( atoi(argv[6])/100., atoi(argv[7])/100., atoi(argv[8])/100. );
  //wf->set_three( atoi(argv[6])/100., atoi(argv[9])/100., atoi(argv[10])/100., atoi(argv[7])/100., atoi(argv[8])/100. );
  //wf->set_four( atoi(argv[6])/100. );
  //wf->set_xxyz( 0 );

  //wf->set_td( atoi(argv[9])/100., atoi(argv[10])/100. );
  //wf->set_r2( atoi(argv[7])/100. );

  //wf->print();

  //create a vmc object and assign the wf
  vmc* myvmc = new vmc();
  myvmc->set_wf( wf );

  //myvmc.thermalize(); //this is done in run()

  //wf->printalpha();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";
  //wf->getwf();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";

  //cout << "accepted: " << myvmc.accepted() << "\n";

//  for(double tt=-1.5; tt<=1.51; tt+=.1)
//  {
    double *mm = new double[NS];
    mm[0]=0.; mm[1]=0.; mm[2]=0.;

    wf->set_hopping( atoi(argv[3])/100., atoi(argv[4])/100., mm );

    wf->findmu( atoi(argv[2]) );

    //wf->print();

    wf->create_ad();

    //wf->set_random_conf();
    //wf->load_alpha();
    //wf->save_alpha();
    //wf->print_alpha();
    //wf->getwf();

    myvmc->initialize( 10 ); //number of bins to average over

    //myvmc.thermalize();

    myvmc->run();

    myvmc->calculate_statistics();

    wf->insert_db();
//  }

  delete[] mm;
  delete myvmc;
  delete wf;

  return 0;
}

