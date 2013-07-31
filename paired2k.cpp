#include "paired2k.h"

paired2k::paired2k(int l, int q) : twoflavor( l, q )
{
  //hopping amplitude (k-space, translational symmetry preserved)
  xi = new double[N];

  //pairing amplitude (k-space, real or complex symmetric for singlet, anti-symmetric for triplet pairing)
#if WFC
  delta = createcomplex( N );
#else
  delta = createdouble( N );
#endif

  N0 = new double[NS];
}

paired2k::~paired2k()
{
  cout << "paired2k::~paired2k()\n";
  delete[] xi;

  destroy(delta, N);
  delete[] N0;
}

//initialize the pairing matrix
void paired2k::create_ad()
{
  //create_dd();

  double xim, en, nx;
  double apx, apy, r;
  int ix,iy,jx,jy,k1x,k1y,k2x,k2y;
  int paires=0;

#if WFC
  complex<double> fk;
#else
  double fk;
#endif

  cout << "paired2k::create_ad() with apx = " << bpars->apx << endl;

  //cout << "delta:\n";
  //write_m(delta, N);

  if(bpars->apx) apx = .5;
  else apx = .0;

  if(bpars->apy) apy = .5;
  else apy = .0;

  for(int j=0; j<N*N; j++) adx[j] = 0.;

  nx = 0.;
  for(int k1=0; k1<N; k1++)
  {
    k1x = (k1%L); k1y = (k1/L);

    for(int k2=0; k2<N; k2++)
    {
      if( abs(delta[k1][k2])<5e-4 ) continue;

      k2x = (k2%L); k2y = (k2/L);

      xim = (xi[k1] + xi[k2])/2.;
      en = sqrt( pow(abs(delta[k1][k2]), 2.) + pow(xim, 2.) );

      fk = delta[k1][k2] / ( en + xim );

      for(int i=0; i<N; i++)
      {
        ix = (i%L); iy = (i/L);

        for(int j=0; j<N; j++)
        {
          jx = (j%L); jy = (j/L);

          r = (double)(ix*k1x + iy*k1y + jx*k2x + jy*k2y) + apx*(ix+jx) + apy*(iy+jy);
#if WFC
          adx[i*N+j] += fk * exp( tPiL*I*r );
#else
          adx[i*N+j] += fk * cos( tPiL*r );
          //adx[i*N+j] += fk * sin( tPiL*r );
#endif
          if( (i==0) && (j==0) ) {
            nx += abs( (1. - xim/en)/2. );
            paires++;
          }
          //paires++;
        }
      }
      //if( paires!= N ) cout << "paires: " << paires << "\n";
   }
  }
  N0[0] = nx; N0[1] = nx;

  if(paires!=N) {
    cout << "WARN: "<< paires << " paired sites...\n";
    exit(-1);
  }
  cff = 1.;
  //normalize
//  for(int i=0; i<N; i++)
//    for(int j=0; j<N; j++) adx[i][j] *= cff/(double)L;

  //cout << "adx:\n";
  //write_m(adx, N);

  cout << "Avg flavor number before projection: ";
  for(int n=0; n<NS; n++) cout << N0[n] << ", ";
  cout << "\n";
}

