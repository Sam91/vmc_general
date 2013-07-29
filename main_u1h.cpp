#include <iostream>
#include <iomanip>

#include "mysql_wrapper.h"

#include "u1.h"
//#include "u1hybrid.h"

#include "vmc.h"


int main(int argc, char *argv[])
{
  int L = atoi(argv[1]);
  
  //u1hybrid* wf = new u1hybrid( L );
  u1* wf = new u1( L*L, L );

  //some preparatory configuration of the wave function

  wf->set_lattice( "triangular" );

  //create a vmc object and assign the wf
  vmc myvmc;
  myvmc.set_wf( wf );

  //myvmc.thermalize(); //this is done in run()

  //wf->printalpha();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";
  //wf->getwf();
  //cout << "wf = " << std::scientific << setprecision(2) << wf->wf << "\n";

  //cout << "accepted: " << myvmc.accepted() << "\n";

//  for(double tt=-1.5; tt<=1.51; tt+=.1)
//  {
    wf->set_hopping( 1., 0. );
//    wf->set_h( 0. );

    wf->create_ad();

    wf->set_random_conf();
    //wf->load_alpha();
    //wf->save_alpha();
    //wf->print_alpha();
    //wf->getwf();

    myvmc.initialize( 20 ); //number of bins to average over

    //myvmc.thermalize();

    myvmc.run();

    myvmc.calculate_statistics();

    //wf->insert_db();
//  }

  return 0;
}
