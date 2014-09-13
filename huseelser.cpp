#include "huseelser.h"

//abstract Huse-Else wave function

huseelser::huseelser( int l ) : wavefunction( l )
{
  //create the d vectors (defining the product state)
#if WFC
  d = createcomplex(NS, N); // N is the number of sites
#else
  d = createdouble(NS, N);
#endif
}

huseelser::~huseelser()
{
  destroy(d, NS);
}

void huseelser::backup_data()
{
  wf_old = wf;
}

void huseelser::print_d()
{
  cout << "d = {";
  for(int i=0; i<N; i++)
  {
    for(int j=0; j<NS; j++) cout << d[j][i] << " ";
    cout << "; ";
  }
  cout << "}";
}

void huseelser::restore_data()
{
  wf = wf_old;
}

//get the wave function wf = <alpha|psi>
void huseelser::getwf()
{
  //product state wave function is given by a complex d = u + iv  vector per site

  wf = 1.;

  //cout << "wf norm: " << wf << "\n";

  //loop over all sites and calculate the product state
  for(int i1=0; i1<N; i1++)
  {
    wf *= d[ alpha->lconf[i1] ][i1];
  }

/*
#if WFC
  wf *= (jastrow3()*dogleg()*cff);
#else
  wf *= (jastrow3()*cff);
#endif
*/
  wf *= (jastrow()*cff);

  //cout << "huseelser::getwf(): " << std::scientific << wf << "\n";

  if( wf != wf )
  {
    cout << "Wf is nan\n";
    throw 4;
  }
}

//swap states on two sites i1 and i2 (corresponding to virtual_replacement)
#if WFC
complex<double> huseelser::swap(int i1, int i2, bool ratio)
{
  complex<double> r;
  
#else
double huseelser::swap(int i1, int i2, bool ratio)
{
  double r;
  
#endif

  if( alpha->lconf[i1] == alpha->lconf[i2] )
  {
    if( !ratio ) cout << "identical states in swap\n";
    r = 1.;
  } else
  {
    if( abs( d[ alpha->lconf[i1] ][i1])>SMALL && abs(d[ alpha->lconf[i2] ][i2])>SMALL )
    {
/*
#if WFC
      r = (d[ alpha->lconf[i2] ][i1]*d[ alpha->lconf[i1] ][i2])/(d[ alpha->lconf[i1] ][i1]*d[ alpha->lconf[i2] ][i2]) * jastrow3(i1, i2) * dogleg(i1, i2);
#else
      r = (d[ alpha->lconf[i2] ][i1]*d[ alpha->lconf[i1] ][i2])/(d[ alpha->lconf[i1] ][i1]*d[ alpha->lconf[i2] ][i2]) * jastrow3(i1, i2);
#endif //WFC
*/
      r = (d[ alpha->lconf[i2] ][i1]*d[ alpha->lconf[i1] ][i2])/(d[ alpha->lconf[i1] ][i1]*d[ alpha->lconf[i2] ][i2]) * jastrow(i1, i2);
      //r = (d[ alpha->lconf[i2] ][i1]*d[ alpha->lconf[i1] ][i2])/(d[ alpha->lconf[i1] ][i1]*d[ alpha->lconf[i2] ][i2]);

    } else {
      cout << "Hitting node swap\n";
      r = 1.;
    }
  }
  
  if( !ratio )
  {
    wf *= r;
    int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
  }
    
  //if( abs(r-1.) > .01 )
  //  cout << "swap (" << alpha->lconf[i1] << ", " << alpha->lconf[i2] << ") returning " << r << "\n";
  return r;
}

#if WFC
complex<double> huseelser::crop(int i1, int i2, bool ratio) {
  complex<double> r;
#else
double huseelser::crop(int i1, int i2, bool ratio) {
  double r;
#endif

  if( alpha->lconf[i1]==0 && alpha->lconf[i2]==0 ) // 0,0 -> 1,2
  {
    if( abs(d[0][i1])>SMALL && abs(d[0][i2])>SMALL ) {
      r = (d[1][i1]*d[2][i2])/(d[0][i1]*d[0][i2]) * jastrow(i1,i2);
    } else {
      //cout << "Hitting node crop\n";
      r = 1.;
    }

  } else { // 1,2 -> 0,0; 2,1 -> 0,0

    if( abs(d[alpha->lconf[i1]][i1])>SMALL && abs(d[alpha->lconf[i1]][i2])>SMALL ) {

      r = (d[0][i1]*d[0][i2])/(d[ alpha->lconf[i1] ][i1]*d[ alpha->lconf[i1] ][i2]) * jastrow(i1,i2);

    } else {
      cout << "Hitting node crop\n";
      r = 1.;
    }
  }

  //cout << "crop returning " << r << "\n";
  return r;
}

void huseelser::find_max_conf() {}

//correct the wf normalization by pow(10,s 20), depending on the sign of s
void huseelser::correct_cff(bool s)
{
  double r;

  if( s )
    r = pow(1.02,N);
  else
    r = pow(.98,N);

  normalize( r );
  //cout << "TF: cff changed to " << cff[n_ext] << endl;
}

void huseelser::normalize(double r)
{
  cff *= r;
  if( std::isinf(cff) )
  {
    cout << "Overflow in cff. resetting...";
    r = 1.;
    cff = 1.;
  }
  cout << "Normalizing with r = " << std::scientific << setprecision(2) << r << "; cff = " << cff << endl;
}

