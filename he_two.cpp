#include "he_two.h"

//define a Huse-Elser wave function for two flavors

he_two::he_two(int n) : huseelser( n )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N; pars->phi1 = 0.; pars->phi2 = 0.;
}

he_two::~he_two()
{
  bpars = NULL;
  delete pars;
}

void he_two::print()
{
  cout << "wf pars->\n";
  cout << " j1, j2 = " << js[0] << ", " << js[1] << "\n";
  cout << " d[n][i] =\n";
  for(int i=0; i<N; i++)
  {
    cout << "   [" << d[0][i] << ", " << d[1][i] << ", " << d[2][i] << "]\n";
  }
}

//set a spiral state; the angles are in units of 2Pi (two dimensions)
void he_two::set_spiral2(double phi1, double phi2)
{
  cout << "he_two::set_spiral2(" << phi1 << ", " << phi2 << ")\n";

  //this only works without sublattice
  if( Q>1 )
  {
    cout << "ERROR: cannot use this function if Q>1\n";
    exit(-1);
  }

  pars->phi1 = phi1;
  pars->phi2 = phi2;
  pars->desc = "he-spiral2";

  int i1=0;
  double cp, sp;

  for(int i=0; i<L; i++)
  {
    for(int j=0; j<L; j++)
    {
      cp = cos( (i*phi1 + j*phi2)*tPi/2. ); //divide by 2 for spin-1/2
      sp = sin( (i*phi1 + j*phi2)*tPi/2. );
#if WFC
      d[0][i1] = complex<double>(cp,sp); //x-y plane
      d[1][i1] = complex<double>(cp,-sp);
#else
      d[0][i1] = cp + sp; //x-z plane
      d[1][i1] = -cp + sp;
#endif
      i1++;
    }
  }
}

//set a spiral state; the angles are in units of 2Pi (one dimensions)
void he_two::set_spiral1( double phi )
{
  cout << "he_two::set_spiral1(" << phi << ")\n";

  //this only works without sublattice
  if( Q>1 )
  {
    cout << "ERROR: cannot use this function if Q>1\n";
    exit(-1);
  }

  ostringstream out;
#if WFC
  out << "he-spiral1c "<< setprecision(1) << phi;
#else
  out << "he-spiral1r "<< setprecision(1) << phi;
#endif

  pars->phi1 = phi;
  pars->desc = out.str();

  int i1=0;
  double cp, sp;

  for(int i=0; i<L; i++)
  {
    cp = cos( ( i*phi )*tPi/2. ); //divide by 2 for spin-1/2
    sp = sin( ( i*phi )*tPi/2. );
#if WFC
    d[0][i1] = complex<double>(cp,sp); //x-y plane
    d[1][i1] = complex<double>(cp,-sp);
#else
    d[0][i1] = cp + sp; //x-z plane
    d[1][i1] = -cp + sp;
#endif
    i1++;
  }
}

//set a four sublattice distortion of the collinear state
void he_two::set_four(double phi1, double phi2)
#if WFC
{
  //this only works without sublattice
  if( Q>1 )
  {
    cout << "ERROR: cannot use this function if Q>1\n";
    exit(-1);
  }

  if( (L%2)!=0 )
  {
    cout << "ERROR: set_four() - even number of sites needed\n";
    exit(-1);
  } 

  cout << "he_two::set_four(" << phi1 << ", " << phi2 << ")\n";
  pars->phi1 = phi1;
  pars->phi2 = phi2;
  pars->desc = "he-four";

  int k, i1=0;
  double cp1, sp1, cp2, sp2;
  cp1 = cos(tPi*phi1/2.);
  sp1 = sin(tPi*phi1/2.);
  cp2 = cos(tPi*phi2/2.);
  sp2 = sin(tPi*phi2/2.);
  
  for(int i=0; i<L; i++)
  {
    for(int j=0; j<L; j++)
    {
      k = (i%2)+2*(j%2);
      switch( k )
      {
        case 0: // A
          d[0][i1] = complex<double>(cp1,sp1); //up (rotated x)
          d[1][i1] = complex<double>(cp1,-sp1); //dn
          break;
        case 1: // B
          d[0][i1] = complex<double>(cp1,-sp1);
          d[1][i1] = complex<double>(cp1,sp1);
          break;
        case 2: // D
          d[0][i1] = complex<double>(cp1+sp1*sp2,sp1*cp2); //up (rotated y)
          d[1][i1] = -complex<double>(cp1-sp1*sp2,-sp1*cp2); //dn
          break;
        case 3: // C
          d[0][i1] = complex<double>(cp1+sp1*sp2,-sp1*cp2);
          d[1][i1] = -complex<double>(cp1-sp1*sp2,sp1*cp2);
          break;
        default :
          cout << "Error\n";
          break;
      }
      i1++;
    }
  }
}
#else
{
  cout << "ERROR: set_four() - complex wf needed\n";
  exit(-1);
}
#endif

//Set the q=0 state on the kagome lattice
void he_two::set_q0()
{
  if( Q!=3 || DIM!=2 )
  {
    cout << "ERROR: we need Q=3 and a Kagome lattice to use he_two::set_q0()\n";
    exit(-1);
  }

  cout << "he_two::set_q0()\n";

#if WFC
  pars->desc = "he-q0";
#else
  pars->desc = "he-q0-r";
#endif

  double cp, sp;

  for(int j=0; j<N; j++)
  {
    switch( j%3 )
    {
      case 0: //A
        cp = 1.;
        sp = 0.;
        break;
      case 1: //B
        cp = cos(tPi/6.);
        sp = sin(tPi/6.);
        break;
      case 2: //C
        cp = cos(tPi/3.);
        sp = sin(tPi/3.);
        break;
      default:
        cp = 0;; sp = 0.;
    }

#if WFC
    d[0][j] = complex<double>(cp,sp); //x-y plane
    d[1][j] = complex<double>(cp,-sp);
#else
    d[0][j] = cp + sp; //x-z plane
    d[1][j] = -cp + sp;
#endif
  }
}

//Set cuboc-1 state on kagome lattice
void he_two::set_cbc1()
{
  if( Q!=3 || DIM!=2 || L%2==1 )
  {
    cout << "ERROR: we need Q=3 and a Kagome lattice for he_two::set_cbc()\n";
    exit(-1);
  }
  cout << "he_two::set_cbc1()\n";

#if WFC
  pars->desc = "he-cbc1";
#else
  cout << "WARN: cannot set cuboc for real wf!\n";
  pars->desc = "he-cbc1-r";
#endif

  int q, j1, j2;
  double phi, theta;
  double v1=0.;
  double v2=0.;
  double v3=0.;

  for(int j=0; j<N; j++)
  {
    q = j % Q; j1 = (j/Q)%L; j2 = (j/Q)/L;

    switch( q )
    {
      case 0:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = -1.; v2 = 1.; v3 = 0.; //12
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = 1.; v2 = -1.; v3 = 0.; //9
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = -1.; v2 = -1.; v3 = 0.; //4
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = 1.; v2 = 1.; v3 = 0.; //1
        }
        break;

      case 1:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = 1.; v2 = 0.; v3 = -1.; //5
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = -1.; v2 = 0.; v3 = -1.; //7
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = 1.; v2 = 0.; v3 = 1.; //10
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = -1.; v2 = 0.; v3 = 1.; //2
        }
        break;

      case 2:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = 0.; v2 = 1.; v3 = 1.; //11
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = 0.; v2 = -1.; v3 = 1.; //6
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = 0.; v2 = -1.; v3 = -1.; //8
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = 0.; v2 = 1.; v3 = -1.; //3
        }
        break;

      default:
        v1 = 0.; v2 = 0.; v3 = 0.;
        cout << "WARN: set_cbc1() error\n";
    }

    phi = atan2( v2, v1 );
    theta = acos( v3/sqrt(v1*v1+v2*v2+v3*v3) );

    //cout << "(" << j1 << ", " << j2 << "); " << q << ": ";
    //cout << "phi= "<< phi/M_PI << "; theta= " << theta/M_PI << ". S = ("<< sin(theta)*cos(phi) << ", " << sin(theta)*sin(phi) << ", " << cos(theta) << ")\n";

#if WFC
    d[0][j] = complex<double>(cos(phi/2.),  sin(phi/2.))*cos(theta/2.);
    d[1][j] = complex<double>(cos(phi/2.), -sin(phi/2.))*sin(theta/2.);
#else
    d[0][j] = cos(phi/2.)*cos(theta/2.);
    d[0][j] = cos(phi/2.)*sin(theta/2.);
#endif
  }
}

//Set cuboc-2 state on kagome lattice
void he_two::set_cbc2()
{
  if( Q!=3 || DIM!=2 || L%2==1 )
  {
    cout << "ERROR: we need Q=3 and a Kagome lattice for he_two::set_cbc2()\n";
    exit(-1);
  }
  cout << "he_two::set_cbc2()\n";

#if WFC
  pars->desc = "he-cbc2";
#else
  pars->desc = "he-cbc2-r";
  cout << "WARN: cannot set cuboc for real wf!\n";
#endif

  int q, j1, j2;
  double phi, theta;
  double v1=0.;
  double v2=0.;
  double v3=0.;

  for(int j=0; j<N; j++)
  {
    q = j % Q; j1 = (j/Q)%L; j2 = (j/Q)/L;

    switch( q )
    {
      case 0:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = -1.; v2 = 1.; v3 = 0.; //12
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = 1.; v2 = -1.; v3 = 0.; //9
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = -1.; v2 = -1.; v3 = 0.; //4
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = 1.; v2 = 1.; v3 = 0.; //1
        }
        break;

      case 1:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = -1.; v2 = 0.; v3 = -1.; //5
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = 1.; v2 = 0.; v3 = -1.; //7
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = -1.; v2 = 0.; v3 = 1.; //10
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = 1.; v2 = 0.; v3 = 1.; //2
        }
        break;

      case 2:
        if( j1%2 == 0 && j2%2 == 0)
        {
          v1 = 0.; v2 = -1.; v3 = -1.; //11
        }

        if( j1%2 == 1 && j2%2 == 0)
        {
          v1 = 0.; v2 = 1.; v3 = -1.; //6
        }

        if( j1%2 == 0 && j2%2 == 1)
        {
          v1 = 0.; v2 = 1.; v3 = 1.; //8
        }

        if( j1%2 == 1 && j2%2 == 1)
        {
          v1 = 0.; v2 = -1.; v3 = 1.; //3
        }
        break;

      default:
        v1 = 0.; v2 = 0.; v3 = 0.;
        cout << "WARN: set_cbc2() error\n";
    }

    phi = atan2( v2, v1 );
    theta = acos( v3/sqrt(v1*v1+v2*v2+v3*v3) );

    //cout << "(" << j1 << ", " << j2 << "); " << q << ": ";
    //cout << "phi= "<< phi/M_PI << "; theta= " << theta/M_PI << ". S = ("<< sin(theta)*cos(phi) << ", " << sin(theta)*sin(phi) << ", " << cos(theta) << ")\n";

#if WFC
    d[0][j] = complex<double>(cos(phi/2.),  sin(phi/2.))*cos(theta/2.);
    d[1][j] = complex<double>(cos(phi/2.), -sin(phi/2.))*sin(theta/2.);
#else
    d[0][j] = cos(phi/2.)*cos(theta/2.);
    d[0][j] = cos(phi/2.)*sin(theta/2.);
#endif
  }
}

// Set the sq3 x sq3 state on the kagome lattice
void he_two::set_sq3()
{
  if( Q!=3 || DIM!=2 )
  {
    cout << "ERROR: we need Q=3 and a Kagome lattice to use he_two::set_sq3()\n";
    exit(-1);
  }

  if( L%3 != 0 )
  {
    cout << "WARN: for sq-3 state, the linear system size should be a multiple of 3 !!!\n\n";
  }

  cout << "he_two::set_sq3()\n";

#if WFC
  pars->desc = "he-sq3";
#else
  pars->desc = "he-sq3-r";
#endif

  int q, j1, j2;
  double cp, sp;

  for(int j=0; j<N; j++)
  {
    q = j % Q; j1 = (j/Q)%L; j2 = (j/Q)/L;

    switch( (j1+j2)%3 )
    {
      case 0: //A
        if( q==1 )
        {
          cp = 1.; sp = 0.;
        } else 
        {
          cp = cos(tPi/6.); sp = sin(tPi/6.);
        }
        break;

      case 1: //B
        if( q==1 )
        {
          cp = cos(tPi/6.); sp = sin(tPi/6.);
        } else 
        {
          cp = cos(tPi/3.); sp = sin(tPi/3.);
        }
        break;

      case 2: //C
        if( q==1 )
        {
          cp = cos(tPi/3.); sp = sin(tPi/3.);
        } else 
        {
          cp = 1.; sp = 0.;
        }
        break;

      default:
        exit(-1);
    }

#if WFC
    d[0][j] = complex<double>(cp,sp); //x-y plane
    d[1][j] = complex<double>(cp,-sp);
#else
    d[0][j] = cp + sp; //x-z plane
    d[1][j] = -cp + sp;
#endif
  }
}

int he_two::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO liquid (id, sites, lattice, txtdsc, NS, N0, mc_length, nbin, ";
  os << "xi1, xi2, xi3, ";
  os << "P1, dP1, P2, dP2, P3, dP3 ) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << mc_length << ", " << nk << ", ";
  os << "round(" << js[0] << ",3), round(" << js[1] << ",3), round(" << js[2] << ",3), ";

  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

int he_two::insert_file(const char *name)
{
  double* f = new double[9];

  for(int i=0; i<3; i++) f[i]   = js[i];
  for(int i=0; i<3; i++) f[i+3] = average[i];
  for(int i=0; i<3; i++) f[i+6] = sigma[i];

  int r = fappend(f, 9, name);

  delete[] f;
  return r;
}

/*
int he_two::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, t1, t2, t3, t1b, t2b, t1c, t2c, jl, phi1, phi2, mc_length, nk, ";
  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4, P1p, dP1p) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", ";
  os << "round(" << js[0] << ",3), round(" << js[1] << ",3), round(" << js[2] << ",3), round(" << js[3] << ",3), round(" << js[4] << ",3), round(" << js[5] << ",3), round(" << js[6] << ",3), round(" << jl << ",3), ";
  os << "round(" << pars->phi1 << ",6), round(" << pars->phi2 << ",6), " << mc_length << ", " << nk << ", ";

  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ", " << average[4] << ", " << sigma[4];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}
*/

//overwrite the default jastrow factors for the anisotropic triangular lattice

//accumulate the nearest neighbor Jastrow factor for this state
/*
double he_two::jastrow()
{
  int nJ1a=0;
  int nJ1b=0;
  int nJ1c=0;
  int nJ2a=0;
  int nJ2b=0;
  int nJ2c=0;
  int nJ3=0;
  
  for(int i1=0; i1<N; i1++)
  {
    //accumulate permutation invariant Jastrow factors
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][1]] ) nJ1a += 2; //add two when same flavor occupies
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][2]] ) nJ1b += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][3]] ) nJ1c += 2;

    nJ1a-=1; //subtract 1 per link
    nJ1b-=1;
    nJ1c-=1;

    //next neighbor
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][2]] ) nJ2a += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][1]] ) nJ2b += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][3]] ) nJ2c += 2;

    nJ2a-=1;
    nJ2b-=1;
    nJ2c-=1;

    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)//third-neighbor
    {
      if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][2][n]] ) nJ3 += 2;
      nJ3--;
    }
  }
  
  double r = exp( - (js[0]*nJ1a + js[1]*nJ2a + js[2]*nJ3 + js[3]*nJ1b + js[4]*nJ2b + js[5]*nJ1c + js[6]*nJ2c)/2. );
  if( abs(r)>1e290 )
  {
    cout << "diverging jastrow\n";
    throw 20;
  }
  cout << "JF0 = " << r << "\n";
  return r;
}

double he_two::jastrow(int i1, int i2 )
{
  if( alpha->lconf[i1] == alpha->lconf[i2] ) return 1.;

  int nJ1a=0;
  int nJ1b=0;
  int nJ1c=0;
  int nJ2a=0;
  int nJ2b=0;
  int nJ2c=0;
  int nJ3=0;

  //loop over the nearest neighbors of i1
  if( i2 != alpha->mylattice->connectivity[i1][0][1] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][1]] ) nJ1a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][1]] ) nJ1a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][2] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][2]] ) nJ1b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][2]] ) nJ1b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][3] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][3]] ) nJ1c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][3]] ) nJ1c++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][4] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][4]] ) nJ1a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][4]] ) nJ1a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][5] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][5]] ) nJ1b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][5]] ) nJ1b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][6] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][6]] ) nJ1c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][6]] ) nJ1c++;
  }

  //loop over the nearest neighbors of i2
  if( i1 != alpha->mylattice->connectivity[i2][0][1] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][1]] ) nJ1a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][1]] ) nJ1a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][2] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][2]] ) nJ1b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][2]] ) nJ1b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][3] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][3]] ) nJ1c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][3]] ) nJ1c++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][4] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][4]] ) nJ1a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][4]] ) nJ1a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][5] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][5]] ) nJ1b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][5]] ) nJ1b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][6] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][6]] ) nJ1c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][6]] ) nJ1c++;
  }

  //loop over the next neighbors of i1
  if( i2 != alpha->mylattice->connectivity[i1][1][1] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][1]] ) nJ2a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][1]] ) nJ2a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][2] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][2]] ) nJ2b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][2]] ) nJ2b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][3] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][3]] ) nJ2c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][3]] ) nJ2c++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][4] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][4]] ) nJ2a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][4]] ) nJ2a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][5] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][5]] ) nJ2b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][5]] ) nJ2b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][6] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][6]] ) nJ2c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][6]] ) nJ2c++;
  }

  //loop over the next neighbors of i2
  if( i1 != alpha->mylattice->connectivity[i2][1][1] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][1]] ) nJ2a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][1]] ) nJ2a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][2] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][2]] ) nJ2b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][2]] ) nJ2b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][3] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][3]] ) nJ2c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][3]] ) nJ2c++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][4] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][4]] ) nJ2a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][4]] ) nJ2a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][5] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][5]] ) nJ2b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][5]] ) nJ2b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][6] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][6]] ) nJ2c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][6]] ) nJ2c++;
  }

  //loop over third neighbors of i1
  for(int n=1; n<=alpha->mylattice->connectivity[i1][2][0]; n++)
  {
    if( i2 == alpha->mylattice->connectivity[i1][2][n] ) continue;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3++;
  }

  //loop over third neighbors of i2
  for(int n=1; n<=alpha->mylattice->connectivity[i2][2][0]; n++)
  {
    if( i1 == alpha->mylattice->connectivity[i2][2][n] ) continue;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3++;
  }

  double r = exp( -(js[0]*nJ1a + js[1]*nJ2a + js[2]*nJ3 + js[3]*nJ1b + js[4]*nJ2b + js[5]*nJ1c + js[6]*nJ2c) );
  
  //cout << "JF = " << r << "\n";
  return r;
}
*/

