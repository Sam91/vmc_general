#include "lwave2.h"

lwave2::lwave2(int l) : lwave2(l, 1) {}

lwave2::lwave2(int l, int q) : paired2k( l, q )
{
  pars = new parameters;
  bpars = pars;

  pars->NF = new int[NS];
  for(int n=0; n<NS; n++) pars->NF[n] = N/NS;

  pars->apx = true;
  pars->apy = false;
  pars->N = N;
  pars->phi1 = 0.; pars->phi2 = 0.;
}

lwave2::~lwave2()
{
  delete[] pars->NF;
  delete pars;
}

/*
 * The linear k (or r) is labled as k = ky + L kx; kx = k/L, ky = k%L
 */

void lwave2::create_dd()
{
  double napx, napy;
  if(pars->apx) napx = .5;
  else napx = 0.;

  if(pars->apy) napy = .5;
  else napy = 0.;

  cout << "create_dd with pars->ap = " << pars->apx << ", " << pars->apy << "\n";

  //dispersion on the triangular lattice
  for(int k=0; k<N; k++) xi[k] = -2.*(pars->t1*cos(tPiL*((double)(k%L)+napx)) + pars->t1c*cos(tPiL*((double)(k/L)+napy)) + pars->t1b*cos(tPiL*( (double)((k%L) + (k/L)) + napx+napy)) ) - (pars->t1+pars->t1b+pars->t1c)/3.*pars->mu; //FIXME: add t2

  int mkx, mky, mk1;

  //pairing amplitude
  for(int k1=0; k1<N; k1++)
  {
    mkx = -(int)pars->apx - (k1%L); if( mkx<0 ) mkx += L;
    mky = -(int)pars->apy - (k1/L); if( mky<0 ) mky += L;

    mk1 =  mkx + L*mky;
    if(mk1<0 || mk1>=N) cout << "mk1 = " << mk1 << "\n";

    for(int k2=0; k2<N; k2++)
    {
      if( k2==mk1 )
      {
#if WFC
//      complex<double> dd0 = pars->dd*( sin(tPiL*(k1%L+nap)) + exp(M_PI/3.*I )*sin(tPiL*(k1/L)) + exp(2.*M_PI/3.*I )*sin(tPiL*(k1%L+nap + (k1/L))) ); //pip
      complex<double> dd0 = pars->dd*( cos(tPiL*((double)(k1%L)+napx)) + exp(-2.*M_PI/3.*I )*cos(tPiL*((double)(k1/L)+napy)) + exp(2.*M_PI/3.*I )*cos(tPiL*( (double)(k1%L + (k1/L)) + napx + napy)) ); //did
#else
      double dd0 = pars->dd*( cos(tPiL*(k1%L+napx)) + (1. - pars->phi1)*cos(tPiL*(k1/L+napy)) + (1. - pars->phi2)*cos(tPiL*(k1%L + (k1/L) + napx + napy)) ); //extended anisotropic s-wave
//      double dd0 = pars->dd*( -.5*cos(tPiL*( (double)(k1%L) + napx)) + (1. - pars->phi1)*cos(tPiL*((double)(k1/L) + napy)) -.5*cos(tPiL*((double)(k1%L + (k1/L) + napx+napy))) ); //d-wave-2 on the triangular lattice
//      double dd0 = pars->dd*( -.5*cos(tPiL*(k1%L+nap)) -.5*cos(tPiL*(k1/L + nap)) + (1. - pars->phi1)*cos(tPiL*(k1%L + (k1/L)+2.*nap)) ); //d-wave-3 on the triangular lattice
//      double dd0 = pars->dd*( sin(tPiL*(k1%L+nap)) + .5*sin(tPiL*(k1/L)) - .5*sin(tPiL*(k1%L+nap + (k1/L))) ); //p-wave
//      double dd0 = pars->dd*( sin(tPiL*(k1%L+nap)) - sin(tPiL*(k1/L)) + sin(tPiL*(k1%L+nap + (k1/L))) ); //f-wave
//      double dd0 = pars->dd*( cos(tPiL*(k1%L+nap)) + cos(tPiL*(k1/L)) + pars->phi1*cos(tPiL*(k1%L+nap + (k1/L))) ); //s-wave on the square lattice
//      double dd0 = pars->dd*( cos(tPiL*(k1%L+nap)) - cos(tPiL*(k1/L)) + pars->phi1*cos(tPiL*(k1%L+nap + (k1/L))) ); //d-wave on the square lattice
#endif
        if( abs(dd0)<1e-3 ) {cout << "lwave2: WARN - vanishing pairing: " << k1 << ", " << k2 << "\n"; delta[k1][k2] = 1e-3;}
        else delta[k1][k2] = dd0;
      } else
        delta[k1][k2] = 0.;
    }
  }
  cff = 1.;
//  pars->desc = "d-wave-tr2";
  //pars->desc = "d-wave-tr";
  //pars->desc = "p-wave";
  //pars->desc = "d-wave";
  //pars->desc = "pip";
  //pars->desc = "did2";
  pars->desc = "s-wave";
}

void lwave2::print()
{
  cout << "wf pars->\n";

  cout << "t1  = " << pars->t1  << "\n";
  cout << "t2  = " << pars->t2  << "\n";
  cout << "dd  = " << pars->dd  << "\n";
  cout << "lth = " << pars->lth << "\n";
  cout << "lr1 = " << pars->lr1 << "\n";
  cout << "lr2 = " << pars->lr2 << "\n";
}

int lwave2::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, apx, apy, "; //general values, always present

  os << "mu, t1, t1b, t1c, t2, dd, lth, lr1, lr2, phi1, phi2, mc_length, nk, "; //wave function specific parameters

  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4, P1p, dP1p) VALUES ("; //measured quantities

  os << "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << pars->NF[0] << ", " << (int)pars->apx << ", " << (int)pars->apy << ", ";

  os << "round(" << pars->mu  << ",3), round(" << pars->t1 << ",3), round(" << pars->t1b << ",3), round(" << pars->t1c << ",3), round(" << pars->t2  << ",3), round(" << pars->dd << ",3),";
  os << "round(" << pars->lth << ",3), round(" << pars->lr1 << ",3), round(" << pars->lr2  << ",3), round(" << pars->phi1 << ",6), round(" << pars->phi2<< ", 6), ";
  os << mc_length << ", " << nk << ", ";
 
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ", " << average[4] << ", " << sigma[4];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

//this function tries to find an apropriate chemical potential mu[0] such that the avg flavor number before projections has N[0]=N[1]=N/2
void lwave2::findmu() { findmu(N2); }
/*{
  int n, nmax;
  cout << "Entering findmu();\n";

  double dm = .1; //delta mu

  nmax = 10;
  for(n=0; n<10; n++)
  {
    create_ad();

    if( abs(N0[0]-N2)<.1 ) {
      cout << "Appropriate muz is " << pars->mu << "\n";
      break;
    } else {
      if( N0[0]>N2 ) pars->mu -= dm;
      else pars->mu += dm;
    }
  }

  if(n==nmax) cout << "No apropriate mu_z found...\n";
}*/

void lwave2::findmu( int nz ) //tries to find a mu[0] such that Nz = nz before projection
{
  int n, nmax;
  cout << "Entering findmu(" << nz << ");\n";

  pars->NF[0] = nz;
  for(int i=1; i<NS; i++) pars->NF[i] = (L2-nz)/(NS-1);

  double dm = .1; //delta mu

  double delta_n, last_mu2, last_mu, last_d2, last_d;
  last_mu = pars->mu-1.;

  nmax = 30;
  for(n=0; n<nmax; n++)
  {
    cout << " muz = " << pars->mu << "\n";
    create_dd();
    create_ad();

    delta_n = abs(N0[0]-nz);
    if( delta_n<.1 )
    {
      cout << "Appropriate muz is " << pars->mu << "\n";
      break;
    } else
    {
      last_mu2 = last_mu; last_mu = pars->mu; //save the last muz
      last_d2 = last_mu2; last_d = delta_n;

      if( N0[0]>nz ) pars->mu -= dm;
      else pars->mu += dm;

      if( abs(last_mu2-pars->mu)<1e-10 )
      {
        cout << "Switching between two mus: ";
        if(last_d2 < last_d)
          pars->mu = last_mu2;
        else
          pars->mu = last_mu;
        cout << "Let's take muz = " << pars->mu << "\n";
        break;
      }
    }
  }

  if(n==nmax) cout << "No apropriate muz found...\n";
}

