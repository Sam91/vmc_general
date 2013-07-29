#include "amperean.h"


amperean::amperean(int n, int l) : paired2k( n, l )
{
  pars = new parameters;
  bpars = pars;

  pars->NF = new int[NS];
  for(int n=0; n<NS; n++) pars->NF[n] = N/NS;

  pars->apx = false;
  pars->apy = false;
  pars->N = N;
  pars->r = 0;
  eps = 1e-3;

  pars->t1=1.;
  pars->t2=0.;
  pars->dd=1.;
  pars->dd0=.01;
  pars->lth=0.;
  pars->lr1=0.;
  pars->lr2=0.;
}

amperean::~amperean()
{
  cout << "amperean::~amperean()\n";
  delete[] pars->NF;
  bpars = NULL;
  delete pars;
}

/*
 * The linear k (or r) is labled as k = ky + L kx; kx = k/L, ky = k%L
 */
/*
void amperean::create_dd()
{
  if( !(pars->r==0 || pars->r==1) ) {
    cout << "ERROR: invalid parameter r\n";
    exit(-1);
  }
  cout << "eps = " << eps << "\n";

  //dispersion on the triangular lattice
  for(int k=0; k<N; k++) xi[k] = - 2.*pars->t1*(cos(tPiL*(k%L+(double)pars->ap/2.)) + cos(tPiL*(k/L)) + cos(tPiL*(k%L+(double)pars->ap/2. + (k/L))) ) - pars->mu; //FIXME: add t2

  int mkx, mky, mk1;
  double k, kf, th0, th, dth, k2;

  //specify the pairing amplitude

  for(int k1=0; k1<N; k1++)
    for(int k2=0; k2<N; k2++)
      delta[k1][k2] = 0.;

  for(int k1=0; k1<N; k1++)
  {
    get_fflo(mkx, mky, k1%L, k1/L, pars->r); //get the fflo paired vector
    get_polar(k, th, 2*(k1%L)+(int)pars->ap, 2*(k1/L)); //get the polar coordinates of k and its partner
    get_polar(k2, th2, 2*mkx+(int)pars->ap, 2*mky);

    if( abs(k-k2)>1e-2 ) cout << "ERROR: mismatch in norm: " << k << ", " << k2 << "\n";
    if( abs(th-th2)>M_PI/3+eps && abs(th-th2)<5./6.*M_PI ) cout << "ERROR: mismatch in angle: " << th << ", " << th2 << "(" << k1%L << ", " << k1/L << "); " << mkx << ", " << mky << ")\n";

    kf = bisect( &xif, 0., M_PI, th, pars->mu); //get the Fermi vector at this angle
    //cout << "kf = " << kf << " for th = " << th << "\n";

    dth = abs( (th-th2)/2. );
    if( dth>M_PI/6.+eps ) dth = M_PI - dth;
    if( dth<-eps || dth>M_PI/6.+eps ) cout << "ERROR: incorrect dy: " << dth << "; (" << k1%L << ", " << k1/L << "); (" << mkx << ", " << mky << "); " << th << ", " << th2 << "\n";

    mk1 =  mkx + L*mky;

    if(k<kf)
      delta[k1][mk1] = pars->dd*exp(-pars->lr1*abs(kf-k)- pars->lth*dth);
    else
      delta[k1][mk1] = pars->dd*exp(-pars->lr2*abs(kf-k)- pars->lth*dth);
  }
  if( pars->r == 0 )
    pars->desc = "amperean-0";
  else
    pars->desc = "amperean-1";

  if( abs(pars->dd) > 1.) cff = 1e-1*pars->dd;
  else cff = 1e-1;
  cff /= (1. + 2.*(pars->lth + pars->lr1 + pars->lr2));
}
*/

void amperean::create_dd()
{
  cout << "amperean::create_dd() with mu = " << pars->mu << endl;

  if( !(pars->r==0 || pars->r==1) )
  {
    cout << "ERROR: invalid parameter r\n";
    exit(-1);
  }

  double kabs, th, th0, kperp, kpara, kf;

//  double dd0 = abs(pars->dd)/100.;
//  if( dd0<1e-3 ) dd0 = 1e-3;

  //dispersion on the triangular lattice (will be used in paired2k)
  for(int k=0; k<N; k++)
    xi[k] = xik( k%L+(double)pars->apx/2., k/L+(double)pars->apy/2. );

  for(int k1=0; k1<N; k1++)
    for(int k2=0; k2<N; k2++) delta[k1][k2] = 0.;

  int kx, ky, k2x, k2y, q1, q2, dx, dy, nkmax, nkmin;

  double xi_n, xi_min;
  int *nfx = new int[7];
  int *nfy = new int[7];

  nfx[0] = -(int)pars->apx; nfy[0] = -(int)pars->apy;

  if( pars->r==0 ) 
  {
    nkmax = 2*L/3; //long axis
//    pars->desc = "amp-0";
    pars->desc = "amp-0-r";
  } else
  {
    nkmax = L; //short axis
//    pars->desc = "amp-1";
    pars->desc = "amp-1-r";
  }
  
  for(int n=1; n<=6; n++)
  {
    if( pars->r==0 )
    {
      switch( n )
      {
        case 1:
          dx = 1; dy = 1;
          break;
        case 2:
          dx = 2; dy = -1;
          break;
        case 3:
          dx = 1; dy = -2;
          break;
        case 4:
          dx = -1; dy = -1;
          break;
        case 5:
          dx = -2; dy = 1;
          break;
        case 6:
          dx = -1; dy = 2;
          break;
        default:
          dx = dy = 0;
      }
    } else //r == 1
    {
      switch( n )
      {
        case 1:
          dx = 1; dy = 0;
          break;
        case 2:
          dx = 1; dy = -1;
          break;
        case 3:
          dx = 0; dy = -1;
          break;
        case 4:
          dx = -1; dy = 0;
          break;
        case 5:
          dx = -1; dy = 1;
          break;
        case 6:
          dx = 0; dy = 1;
          break;
        default:
          dx = dy = 0;
      }
    }

    xi_min = 1e6; nkmin = -1;

    for(int nk=0; nk<nkmax; nk++)
    {
      kx = dx*nk; ky = dy*nk;
      if( kx<0 ) kx += 2*L; if( kx>2*L ) kx -= 2*L;
      if( ky<0 ) ky += 2*L; if( ky>2*L ) ky -= 2*L;

      xi_n = xik( (double)kx/2., (double)ky/2. );
      if( (abs(xi_n) < xi_min) && xi_n>=0. )  //take the k-point with smallest POSITIVE dispersion
      {
        xi_min = abs(xi_n); nkmin = nk;
      }
    }
    cout << "Fermi nk = " << nkmin << " for hexant " << n << "\n";

    nfx[n] = dx*nkmin - (int)pars->apx; nfy[n] = dy*nkmin - (int)pars->apy;
    if( nfx[n]<0 ) nfx[n] += 2*L; if( nfx[n]>2*L ) nfx[n] -= 2*L;
    if( nfy[n]<0 ) nfy[n] += 2*L; if( nfy[n]>2*L ) nfy[n] -= 2*L;
  }

  //double kf = sqrt(dot(nfx[1], nfx[1])); //suppose that all Q vectors have the same distance to the center of the BZ
  int namp = 0; //count the number of states paired in the amperean channel
  double dkx=0.;
  double dky=0.;
  
  for(int k1=0; k1<N; k1++) //loop over all k-points
  {
    kx = (k1%L); ky = (k1/L);

    q1 = getQ(kx, ky, pars->r); //get the hexant of this k-point

    k2x = nfx[q1] - kx; k2y = nfy[q1] - ky; //get the inversion-symmetric point

    q2 = getQ(k2x, k2y, pars->r);

    if( q1!=q2 ) //if the reverted point is outside, then reflect at the origin
    {
      q1 = 0;
      k2x = nfx[q1] - kx; k2y = nfy[q1] - ky;
    }

    if( k2x < 0 )  k2x += L;
    if( k2y < 0 )  k2y += L;
    if( k2x >= L ) k2x -= L;
    if( k2y >= L ) k2y -= L;
    int k2 = k2x + L*k2y;

    //if xibar<0, reflect at 0 instead (to avoid non-analytic pairing)
//    if( (q1!=0) && (xi[k1]+xi[k2]<=0.) )
    if( (q1!=0) && ((xi[k1]<0.) || (xi[k2]<0.)) )
    {
      q1 = 0;
      k2x = nfx[q1] - kx; k2y = nfy[q1] - ky;
      if( k2x < 0 )  k2x += L;
      if( k2y < 0 )  k2y += L;
      if( k2x >= L ) k2x -= L;
      if( k2y >= L ) k2y -= L;
      k2 = k2x + L*k2y;
    }

    if( q1==0 ) //case of boundary or inner k-point, pair in s-wave
      delta[k1][k2] = pars->dd0;
    else
    {
      if( abs(pars->lth)<1e-3 )
        delta[k1][k2] = pars->dd;
      else
      {
        get_polar(kabs, th, 2*kx+(int)pars->apx, 2*ky+(int)pars->apy);
        get_polar(kf, th0, nfx[q1], nfy[q1]);
        th0 = getQ_angle( q1, pars->r );
      
        kpara = kabs*cos(th-th0);
        kperp = kabs*sin(th-th0);
        //cout << "Found k = (" << (kpara-kf) << ", " << kperp << ")\n";
        if(abs(kpara-kf)>dkx) dkx = abs(kpara-kf);
        if(abs(kperp)>dky) dky = abs(kperp);
        
        if( kpara<0 )
          cout << "ERROR with k = (" << kx << ", " << ky << "), q = " << q1 << ", th = " << th << ", th0 = " << th0 << "\n";

        delta[k1][k2] = pars->dd*pow( pow(kperp,2) + pow(pars->lr1*(kpara-kf),2) + pow(pars->lr2,2), -pars->lth/2. );
      }
      namp++;
    }
  }

  if( dkx>0 || dky>0. )
    cout << "dk-max = (" << dkx << ", " << dky << ")\n";

/*  cff = 5e-2;
  if(L>=16) cff = 1e-2;
  if(L>=26) cff = 8e-3;
  if( abs(pars->dd)>=2. ) cff *= 1.8;
  else if( abs(pars->dd)>1. ) cff *= 1.3;
*/
  cout << "n-amp: " << namp/2 << " ~ " << (double)(100*namp)/(double)(2*L2) << "%\n";

  delete[] nfx; delete[] nfy;
}

/*int amperean::diff(int k1, int k2) //get the smallest difference between two pts in the BZ
{
  int k1x = k1%(2*L);
  int k1y = k1/(2*L);
  int k2x = k2%(2*L);
  int k2y = k2/(2*L);

  k2x -= k1x;
  k2y -= k1y;

  if( k2x >= L ) k2x -= 2*L;
  if( k2y >= L ) k2y -= 2*L;
  if( k2x < L  ) k2x += 2*L;
  if( k2x < L  ) k2y += 2*L;

  return k2x + 2*L*k2y;
}
*/
// Implements the mapping (k1,k2) -> (k,theta) on the triangular lattice
void amperean::get_polar(double &k, double &theta, int kx, int ky)
{
  int k0x, k0y;

//  int kx = 2*k1 + (int)pars->ap;
//  int ky = 2*k2;

  //cout << "dk = (" << dk1 << ", " << dk2 << ")\n";

  //determine the closest reciprocal lattice vector
  if( (ky <= 2*(L-kx)) && (2*ky <= 2*L-kx) )
  {
    k0x = 0; k0y = 0;
  } else
    if( (ky > 2*(2*L-kx)) && (2*ky > 4*L-kx) )
    {
      k0x = 2*L; k0y = 2*L;
    }
    else
      if( ky >= kx )
      {
        k0x = 0; k0y = 2*L;
      } else
      {
        k0x = 2*L; k0y = 0;
      }
  
  //cout << "k0 = (" << k0x << ", " << k0y << ")\n";

  double dkx = .5*sqrt(3.)*(double)(kx - k0x);
  double dky = (double)(ky - k0y) + .5*(double)(kx - k0x);

  k = sqrt(dkx*dkx + dky*dky)/(double)(2*L);
  theta = atan2( dky, dkx );

  //cout << "get_polar(ap=" << pars->ap << "): (" << kx << ", " << ky << ") -> " << k << ", " << theta << "\n";
}

// Inverse mapping from the above
void amperean::get_carth(int& k1, int& k2, double k, double th)
{
  if( (th<-M_PI) || (th>M_PI) )//bring this into the interval [-Pi, Pi] 
  {
    cout << "correcting theta\n";
    th = th/M_PI + 1.;
    th = M_PI*(th - floor(th) - 1.);
  }

  double dk1 = k*cos(th); //these are NOT the square coordinates, but the coords in the undistorted BZ
  double dk2 = k*sin(th);
  double k01, k02;

  //determine the quadrant
  if( (th>=M_PI/6.-SMALL) && (th<=M_PI/2.+SMALL) ) { k01 = 0.; k02 = 0.; } //[pi/6, pi/2]
  else
    if( ((th>M_PI/2.) && (th<=7.*M_PI/6.)) || ((th>=-M_PI) && (th<=-5.*M_PI/6.)) ) { k01 = 1.; k02 = 0.; }
    else
      if( (th>-5.*M_PI/6.) && (th<-M_PI/2.) ) { k01 = 1.; k02 = 1.; }
      else
        if( (th>=-M_PI/2.-SMALL) && (th<M_PI/6.) ) { k01 = 0.; k02 = 1.; }
        else {
          cout << "ERROR: no quadrant found\n"; k01 = k02 = 0.;
        }
  //cout << "k0 = (" << k01 << ", " << k02 << ")\n";

  dk1 *= 2./sqrt(3.); //rotate it to the square
  dk2 = dk2 - .5*dk1;

  dk1 += k01;
  dk2 += k02;

  //cout << "dk = (" << dk1 << ", " << dk2 << ")\n";

  k1 = (int)((double)L*dk1 - (double)pars->apx/2. + .5);
  k2 = (int)((double)L*dk2 - (double)pars->apy/2. + .5);

  //cout << k1 << ", " << k2 << "\n";
  if(k1>=L) k1 -= L; if(k1<0) k1 += L;
  if(k2>=L) k2 -= L; if(k2<0) k2 += L;

  //cout << "get_carth(ap=" << pars->ap << "): (" << k << ", " << th << ") -> " << k1 << ", " << k2 << "\n";
}

//returns the symmetric partner in the BZ at (k,th2), in the 6-fold symmetric pairing (r=0: long symmetry axis; r=1: short symmetry axis)
void amperean::get_fflo(int &k1, int &k2, int j1, int j2, int r)
{
  double k, th, th2, thn;

  get_polar(k, th, 2*j1+(int)pars->apx, 2*j2+(int)pars->apy); //get the momentum verctor in polar coordinates

  th2 = -10.;

  switch( r )
  {
    case 0:
      if( abs(th)>5./6.*M_PI+eps ) th2 = -th;
      else {
        for(int n=-2; n<3; n++) {
          thn = ((double)n)*M_PI/3.;
          //if(n==-2) cout << "("<< j1 <<","<< j2<< "): " << abs(th-thn) << ", " << M_PI/6.+eps << ";"<< (abs(th-thn)<=M_PI/6.+eps) << ", " << abs(n)%2 << "\n";
          if( (((abs(n)%2)==0) && abs(th-thn)<=M_PI/6.+eps) || (((abs(n)%2)==1) && abs(th-thn)<M_PI/6.-eps) ) { th2 = 2.*thn - th; }
        }
      }
      break;

    case 1:
      if( th<=-M_PI+eps ) { th2 = 4./3.*M_PI - th; }
      else if( th>=M_PI-eps ) th2 = 4./3.*M_PI - th + 2.*M_PI;
      else {
        for(int n=-2; n<=3; n++) {
          thn = ((double)n-.5)*M_PI/3.;
          if( ((abs(n)%2==0) && (abs(th-thn)<=M_PI/6.+eps) ) || ((abs(n)%2==1) && (abs(th-thn)<M_PI/6.)) ) { th2 = 2.*thn - th; }
        }
      }
      break;

    default:
      cout << "ERROR: invalid r = " << r << "\n";
  }

  if( th2<-5. ) {
    cout << "ERROR: no theta2 found: (" << j1 << ", " << j2 << "); " << th/M_PI << "\n";
  }

  get_carth(k1, k2, k, th2);
  //cout << "th: " << th << " -> " << th2 << "\n";
  //cout << j1 << ", " << j2 << " -> " << k1 << ", " << k2 << "\n";
}

double amperean::getQ_angle( int q, int r )
{
  int m = 2-q;
  if( m<=-3 ) q += 6;

  switch( r )
  {
  case 0:
    return (double)m * M_PI/3.;
  case 1:
    return (double)(2*m-1) * M_PI/6.;
  default:
    cout << "ERROR; incorrect parameter\n";
    return 0.;
  }
}

//get the hexant for (k1,k2); 1, ..., 6 (conter-clockwise labelling, starting with x=y direction
int amperean::getQ(int k1, int k2, int r)
{
  int kx, ky;
  kx = 2*k1+(int)pars->apx; ky = 2*k2+(int)pars->apy;

  if( kx<0 ) kx += 2*L; if( kx>=2*L ) kx -= 2*L;
  if( ky<0 ) ky += 2*L; if( ky>=2*L ) ky -= 2*L;

  switch( r )
  {
  case 0:

    if( ky == 0 ) return 0;
    if( kx == 0 ) return 0;
    if( ky == 2*L - kx ) return 0;

    if( (ky <= 2*(L-kx)) && (2*ky <= 2*L-kx) )
    {
      if( (ky == 2*(L-kx)) || (2*ky == 2*L-kx) ) return 0;
      else return 1;
    }

    if( (ky >= 2*(2*L-kx)) && (2*ky >= 4*L-kx) )
    {
      if( (ky == 2*(2*L-kx)) || (2*ky == 4*L-kx) ) return 0;
      else return 4;
    }

    if( ky == kx ) return 0;

    if( ky > kx )
    {
      if( ky > 2*L-kx ) return 2;
      else return 3;
    }
    else
    {
      if( ky > 2*L-kx ) return 6;
      else return 5;
    }
    cout << "ERROR\n"; exit(-1);

  case 1:

    if( kx == ky ) return 0;
    if( ky == 2*(L-kx) ) return 0;
    if( ky == 2*(2*L-kx) ) return 0;
    if( 2*ky == 2*L-kx ) return 0;
    if( 2*ky == 4*L-kx ) return 0;

    if( ky < kx )
    {
      if( ky < 2*(L-kx) ) return 1;
      if( 2*ky < 2*L-kx ) return 4;
      if( ky < 2*(2*L-kx) ) return 5;
      if( 2*ky < 4*L-kx ) return 6;
      return 3;
    } else {
      if( 2*ky < 2*L-kx ) return 6;
      if( ky < 2*(L-kx) ) return 3;
      if( 2*ky < 4*L-kx ) return 2;
      if( ky < 2*(2*L-kx) ) return 1;
      return 4;
    }

    break;

  default:
    cout << "Incorrect parameter\n";
  }
  return -1;
}

void amperean::print()
{
  cout << "wf pars->\n";

  cout << "t1  = " << pars->t1  << "\n";
  cout << "t2  = " << pars->t2  << "\n";
  cout << "dd  = " << pars->dd  << "\n";
  cout << "lth = " << pars->lth << "\n";
  cout << "lr1 = " << pars->lr1 << "\n";
  cout << "lr2 = " << pars->lr2 << "\n";
}

int amperean::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, apx, apy, "; //general values, always present

  os << "mu, t1, t2, dd, dd0, lth, lr1, lr2, "; //wave function specific parameters

  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4) VALUES ("; //measured quantities

  os << "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << pars->NF[0] << ", " << (int)pars->apx << ", " << (int)pars->apy << ", ";

  os << "round(" << pars->mu  << ",3), round(" << pars->t1 << ",3), round(" << pars->t2  << ",3), round(" << pars->dd << ",3), round(" << pars->dd0 << ",3), ";
  os << "round(" << pars->lth << ",3), round(" << pars->lr1 << ",3), round(" << pars->lr2 << ", 3), ";

  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

//this function tries to find an apropriate chemical potential mu[0] such that the avg flavor number before projections has N[0]=N[1]=N/2
void amperean::findmu() { findmu(N2); }
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

void amperean::findmu( int nz ) //tries to find a mu[0] such that Nz = nz before projection
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
    //cout << " muz = " << pars->mu << "\n";
    create_dd();
    create_ad();

    delta_n = abs(N0[0]-nz); //particle-number difference

    if( delta_n<.1 ) {
      cout << "Appropriate mu is " << pars->mu << "\n";
      break;
    } else {
      last_mu2 = last_mu; last_mu = pars->mu; //save the last muz
      last_d2 = last_mu2; last_d = delta_n;

      if( N0[0]>nz ) pars->mu -= dm;
      else pars->mu += dm;

      if( abs(last_mu2-pars->mu)<1e-10 )
      {
        cout << "Switching between two mus...\n";
        if(last_d2 < last_d) {
          pars->mu = last_mu2;
          create_dd();
          create_ad();
        } else
          pars->mu = last_mu;
        cout << "Let's take mu = " << pars->mu << endl;
        break;
      }
    }
  }

  if(n==nmax) cout << "No appropriate mu found...\n";

  normalize(5e-2/(double)L); //normalize the wf tentatively
}

double amperean::xif(double k, double theta, double mu)
{
  return -2.*(cos(k*cos(theta)) + 2*cos(k*cos(theta)/2.)*cos(k*sin(theta)*pow(3.,.5)/2.) ) - mu;
}

double amperean::xik(double kx, double ky)
{
  return - 2.*pars->t1*(cos( tPiL*kx ) + cos( tPiL*ky ) + cos( tPiL*( kx+ky ) ) ) - pars->mu; //FIXME: add t2
}

