#include "u1kagome.h"
#include <iomanip>

// Define an unpaired real fermionic state on the Kagome lattice with first, second, and diagonal neighbor hopping
// The state may have unit-cell doubling (epsilon) and an additional rotation staggering of the hopping real (tau_R)

u1kagome::u1kagome(int l) : u1hybrid( l )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N;

  pars->ap = new bool[DIM];

//set the default values
  pars->ap[0] = false; pars->ap[1] = false;

  pars->e2 = false;
  pars->TR = false;
  pars->gR = 0;

  int nxi = 3;
  pars->xi = new double[nxi]; for(int i=0; i<nxi; i++) pars->xi[i] = 0.;
  pars->dd = new double[nxi]; for(int i=0; i<nxi; i++) pars->dd[i] = 0.;
  pars->a  = new double[nxi]; for(int i=0; i<nxi; i++) pars->a[i]  = 0.;
  pars->b  = new double[nxi]; for(int i=0; i<nxi; i++) pars->b[i]  = 0.;
  pars->ll = new double[nxi]; for(int i=0; i<3; i++)   pars->ll[i] = 0.;
}

u1kagome::~u1kagome()
{
  delete[] pars->xi; delete[] pars->dd;
  delete[] pars->a; delete[] pars->b;
  delete[] pars->ap;

  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1kagome::set_hopping()
{
  u1hybrid::set_hopping(1., 0, pars->ap[0]);
}

//set a real t1-t2-t3 hopping state on Kagome 
//set the real hopping matrix
void u1kagome::set_hoppingk(double mu0)
{
  for(int n=0; n<NS; n++) mu[n] = mu0;

  int i2;
  int q1, q2;

  cout << "ap = ["<< pars->ap[0] << "; " << pars->ap[1] << "]" << endl;
  cout << "e2 = "<< pars->e2 << endl;
  cout << "TR = "<< pars->TR << endl;
  cout << "gR = "<< pars->gR << endl;

  cout << "xi = ["<< pars->xi[0] << "; " << pars->xi[1] << "; " << pars->xi[2] << "]" << endl;
#if WFC
  cout << " a = ["<< pars->a[0] << "; " << pars->a[1] << "; " << pars->a[2] << "]" << endl;
#endif

  if( pars->gR<0 || pars->gR>3 ) {
    cout << "ERROR: invalid value of gR\n";
    exit(-1);
  }

  bool rr;
  if( pars->gR==0 || pars->gR==3 ) rr = pars->TR;
  else rr = !pars->TR;

  double srr; if( rr ) srr = -1.; else srr = 1.;
#if WFC
  double sri; if( pars->TR ) sri = -1.; else sri = 1.;
#endif


  //set the hopping matrix to zero
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int nn1=0; nn1<NS; nn1++)
        for(int nn2=0; nn2<NS; nn2++) t[nn1][nn2][i][j] = 0.;

  int *n1 = new int[DIM];
  int *n2 = new int[DIM];

  for(int i1=0; i1<N; i1++) //loop over all sites
  {
    alpha->mylattice->getnq(n1, q1, i1); //get the carthesian coordinates for i1.

    for(int n=0; n<NS; n++) //loop over flavors (two)
    {
      //loop over first neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][0][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][0][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); // get the cartesian coordinates for i2

        // set hopping on nearest neighbors according to the rotation sign (note that the chosen link order is important here)
        if( (q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1) ) {
#if WFC
          t[n][n][i1][i2] = pars->xi[0]*complex<double>( srr*cos(M_PI*pars->a[0]), +sri*sin(M_PI*pars->a[0]) );
          t[n][n][i2][i1] = pars->xi[0]*complex<double>( srr*cos(M_PI*pars->a[0]), -sri*sin(M_PI*pars->a[0]) );
        } else {
          t[n][n][i1][i2] = pars->xi[0]*complex<double>( cos(M_PI*pars->a[0]), +sin(M_PI*pars->a[0]) );
          t[n][n][i2][i1] = pars->xi[0]*complex<double>( cos(M_PI*pars->a[0]), -sin(M_PI*pars->a[0]) );
#else
          t[n][n][i1][i2] = t[n][n][i2][i1] = srr*pars->xi[0];
        } else {
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[0];
#endif
        }

        if( pars->e2 ) // unit cell doubling
        {
          if( q1==0 && k==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }

          if( (n1[1])%2==0 ) {
            if( q1==1 && k==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==2 && k==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
        }

        // boundary conditions
        if( pars->ap[0] && ( abs(n1[0] - n2[0]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->ap[1] && ( abs(n1[1] - n2[1]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop first neighbors

      //loop over second neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][1][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][1][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); //get the cartesian coordinates for i2

        //set hopping on second neighbors according to rotation sign
        if( (q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1) ) {
#if WFC
          t[n][n][i1][i2] = pars->xi[1]*complex<double>( srr*cos(M_PI*pars->a[1]), +sri*sin(M_PI*pars->a[1]) );
          t[n][n][i2][i1] = pars->xi[1]*complex<double>( srr*cos(M_PI*pars->a[1]), -sri*sin(M_PI*pars->a[1]) );
        } else {
          t[n][n][i1][i2] = pars->xi[1]*complex<double>( cos(M_PI*pars->a[1]), +sin(M_PI*pars->a[1]) );
          t[n][n][i2][i1] = pars->xi[1]*complex<double>( cos(M_PI*pars->a[1]), -sin(M_PI*pars->a[1]) );
#else
          t[n][n][i1][i2] = t[n][n][i2][i1] = srr*pars->xi[1];
        } else {
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[1];
#endif
        }

        if( pars->e2 ) //unit cell doubling
        {
          if( (n1[1])%2==0 ) {
            if( (q1==0 && k==2) || (q1==1 && k==2) || (q1==2 && k==1) ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==1 && k==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
        }

        // boundary conditions
        if( pars->ap[0] && ( abs(n1[0] - n2[0]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->ap[1] && ( abs(n1[1] - n2[1]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop second neighbors

      //loop over third neighbors (diagonal link)
      for(int k=1; k<=alpha->mylattice->links[i1][2][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][2][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); //get the cartesian coordinates for i2

        //setting hopping to one on all third neighbors
        if( q1!=1 ) {
#if WFC
          t[n][n][i1][i2] = pars->xi[2]*complex<double>( srr*cos(M_PI*pars->a[2]), +sri*sin(M_PI*pars->a[2]) );
          t[n][n][i2][i1] = pars->xi[2]*complex<double>( srr*cos(M_PI*pars->a[2]), -sri*sin(M_PI*pars->a[2]) );
        } else {
          t[n][n][i1][i2] = pars->xi[2]*complex<double>( cos(M_PI*pars->a[2]), +sin(M_PI*pars->a[2]) );
          t[n][n][i2][i1] = pars->xi[2]*complex<double>( cos(M_PI*pars->a[2]), -sin(M_PI*pars->a[2]) );
#else
          t[n][n][i1][i2] = t[n][n][i2][i1] = srr*pars->xi[2];
        } else {
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[2];
#endif
        }

        if( pars->e2 ) //unit cell doubling
        {
          if( (n1[1])%2==0 ) {
            if( q1==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
        }

        // boundary conditions
        if( pars->ap[0] && ( abs(n1[0] - n2[0]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->ap[1] && ( abs(n1[1] - n2[1]) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop third neighbors

    } //end loop flavors
  } //end loop sites

  delete[] n1; delete[] n2;
}

void u1kagome::set_hopping3(double *tt1, double *tt2)
{
  u1hybrid::set_hopping3(tt1, tt2, pars->ap[0]);
}

void u1kagome::print() {}

/*
void u1kagome::print_avgs()
{
  cout<< "ff={";
  for(int no=0; no<NO-1; no++) {
    cout << std::fixed << setprecision(5) << average[no] << ", ";
  }
  cout << std::fixed << setprecision(5) << average[NO-1] << "};"<<endl;
  cout<< "ffsigma={";
  for(int no=0; no<NO-1; no++) {
    cout << std::fixed << setprecision(5) << sigma[no] << ", ";
  }
  cout << std::fixed << setprecision(5) << sigma[NO-1] << "};"<<endl;
}
*/

int u1kagome::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO liquid (id, sites, lattice, txtdsc, NS, N0, apx, apy, mc_length, nbin, e2, TR, gR, ";
  os << "l3, l2, l1, xi1, xi2, xi3, dd1, dd2, dd3, a1, a2, a3, b1, b2, b3, ";
  os << "P1, dP1, P2, dP2, P3, dP3 ) VALUES (";

  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << (int)pars->ap[0] << ", " << (int)pars->ap[1] << ", " << mc_length << ", " << nk << ", ";
  os << (int)pars->e2 << ", " << (int)pars->TR << ", " << (int)pars->gR << ", ";

  os << "round(" << mu[0] << ",3), round(" << pars->ll[2] << ",3), round(" << pars->ll[1] << ",3), ";
  os << "round(" << pars->xi[0] << ",3), round(" << pars->xi[1] << ",3), round(" << pars->xi[2] << ",3), ";
  os << "round(" << pars->dd[0] << ",3), round(" << pars->dd[1] << ",3), round(" << pars->dd[2] << ",3), ";
  os << "round(" << pars->a[0]  << ",6), round(" << pars->a[1]  << ",6), round(" << pars->a[2]  << ",6), ";
  os << "round(" << pars->b[0]  << ",6), round(" << pars->b[1]  << ",6), round(" << pars->b[2]  << ",6), ";

  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

int u1kagome::insert_file(const char *name)
{
  double* f = new double[9];

  for(int i=0; i<3; i++) f[i]   = js[i];
  for(int i=0; i<3; i++) f[i+3] = average[i];
  for(int i=0; i<3; i++) f[i+6] = sigma[i];

  int r = fappend(f, 9, name);

  delete[] f;
  return r;
}

