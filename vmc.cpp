#include "vmc.h"
#include <iomanip>

vmc::vmc()
{
  randomize();
  f = 0;
}

vmc::~vmc()
{
  //if( f != 0 ) destroy(f, nk);
  cout << "destroy with nk=" << nk << "\n"; 
  mywf->destroy_f( nk );
}

void vmc::set_wf( wavefunction* wf)
{
  mywf = wf;
  L = wf->getL(); LD = pow(L,DIM);
  N = wf->getN();

  wf->set_walk_length( 2*N ); //sime integer multiple of N should be OK here (2x or 4x)
  //wf->set_mc_length( 40 );
}

//set the number of bins
void vmc::initialize(int n)
{
  if( n < 2 ) {
    cout << "ERROR: vmc::initialize() bin number must be > 1\n";
    exit(-1);
  }
  nk = n;
  cout << "init with nk=" << nk << "\n"; 
  mywf->initiate_f( nk );
}

void vmc::thermalize()
{
  cout << "Entering thermalization\n";

  mywf->find_starting_conf();

  int therm = 100;

  while( 1 )
  {
    mywf->accepted = 0; cout << "Reset acceptance\n";

    for(int i=0; i<therm; i++)
    {
      if( abs(mywf->wf)>1e250 ) {
        cout << "WARN(vmc): diverging wf=" << mywf->wf << "\n";
        mywf->correct_cff( false ); i=0; mywf->accepted = 0; mywf->getwf();
        //cout << "New wf = " << mywf->wf << endl;
      }
      mywf->walk(); cout << ".";
    }
    cout << " thermalized.\n";

    if( abs(mywf->wf) > 1e-4 ) break;
    else
    {
      cout << "WARN(vmc): vanishing wf=" << mywf->wf << "\n";
      mywf->correct_cff( true );
      mywf->getwf();
      cout << "Rethermalizing \n";
    }
  }

//    cout << "acc " << mywf->step() << "\n";

  cout << "Thermalization acceptance: " << accepted() << "/" << therm*mywf->get_walk_length() << " = ";
  cout << setprecision(1) << std::fixed << 100.*(double)accepted()/(double)(therm*mywf->get_walk_length()) << "%. wf = " << std::scientific << mywf->wf << endl;

  if( accepted()<5 ) exit(-1);
}

int vmc::accepted() { return mywf->accepted; }

void vmc::accumulate()
{
  mywf->accumulate();
}

void vmc::run()
{
  thermalize();

  cout << "Starting run with nk=" << nk << "\n";
  double time_start = time(0);
  mywf->reset_run();

  for(int i=0; i<nk; i++) //loop over bins for variance
  {
    if( i>0 && i%20 == 0 ) estimate_time( time_start, nk, i );

    mywf->walk_accumulate();

    mywf->collect_data();

  }
  estimate_time( time_start, nk, nk);
}

void vmc::calculate_statistics()
{
  mywf->calculate_statistics();
}

