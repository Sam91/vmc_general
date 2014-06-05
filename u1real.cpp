#include "u1real.h"
#include "spinone.h"
#include <iomanip>

// Define an unpaired real fermionic state on the Kagome lattice with first, second, and diagonal neighbor hopping
// The state may have unit-vell doubling (epsilon) and an additional rotation staggering of the hopping real (tau_R)

u1real::u1real(int l, int q) : u1hybrid( l, q )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N;

//set the default values
  pars->apx = false; pars->apy = false;

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

u1real::~u1real()
{
  bpars = NULL;
  delete[] pars->xi; delete[] pars->dd;
  delete[] pars->a; delete[] pars->b;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1real::set_hopping()
{
  u1hybrid::set_hopping(1., 0, pars->apx);
}

//set a real t1-t2-t3 hopping state on Kagome 
//set the real hopping matrix
void u1real::set_hoppingk(double mu0)
{
  for(int n=0; n<NS; n++) mu[n] = mu0;

  int i2;
  int nx2, ny2, q2, nx1, ny1, q1;

  cout << "ap = ["<< pars->apx << "; " << pars->apy << "]" << endl;
  cout << "e2 = "<< pars->e2 << endl;
  cout << "TR = "<< pars->TR << endl;
  cout << "gR = "<< pars->gR << endl;

  cout << "xi = ["<< pars->xi[0] << "; " << pars->xi[1] << "; " << pars->xi[2] << "]" << endl;

  if( pars->gR<0 || pars->gR>3 ) {
    cout << "ERROR: invalid value of gR\n";
    exit(-1);
  }

  bool rr;
  if( pars->gR==0 || pars->gR==3 ) rr = pars->TR;
  else rr = !pars->TR;

  //set the hopping matrix to zero
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n1=0; n1<NS; n1++)
        for(int n2=0; n2<NS; n2++) t[n1][n2][i][j] = 0.;

  for(int i1=0; i1<N; i1++) //loop over all sites
  {
    alpha->mylattice->getnq(nx1, ny1, q1, i1); //get the carthesian coordinates for i1.

    for(int n=0; n<NS; n++) //loop over flavors (two)
    {
      //loop over first neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][0][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][0][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(nx2, ny2, q2, i2); // get the cartesian coordinates for i2

        // set hopping on nearest neighbors according to the rotation sign (note that the chosen link order is important here)
        if( rr && ((q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1)) )
          t[n][n][i1][i2] = t[n][n][i2][i1] = -pars->xi[0];
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[0];

        if( pars->e2 ) // unit cell doubling
        {
          if( q1==0 && k==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }

          if( ny1%2==0 ) {
            if( q1==1 && k==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==2 && k==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
        }

/*
          if( (q1==1) && (q2==0) ) // base link of dn triangle
            if( ny1%2==0 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }

          if( (q1==2) && (q2==0) ) // left link of up triangle
            if( ny1%2==0 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          if( (q1==1) && (q2==2) ) { // right link of up triangle
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }
*/

        // boundary conditions
        if( pars->apx && ( abs(nx1 - nx2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->apy && ( abs(ny1 - ny2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop first neighbors

      //loop over second neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][1][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][1][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(nx2, ny2, q2, i2); //get the cartesian coordinates for i2

        //set hopping on second neighbors according to rotation sign
        if( rr && ((q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1)) )
          t[n][n][i1][i2] = t[n][n][i2][i1] = -pars->xi[1];
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[1];

        if( pars->e2 ) //unit cell doubling
        {
          if( ny1%2==0 ) {
            if( (q1==0 && k==2) || (q1==1 && k==2) || (q1==2 && k==1) ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==1 && k==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
/*
          if( ((q1==1)&&(q2==2))||((q1==2)&&(q2==1)) ) //base link of dn triangle
            if( ny1%2==0 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }

          if( ((q1==1)&&(q2==0))||((q1==2)&&(q2==0)) ) //left link of up triangle
            if( ny1%2==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
*/
        }

        // boundary conditions
        if( pars->apx && ( abs(nx1 - nx2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->apy && ( abs(ny1 - ny2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop second neighbors

      //loop over third neighbors (diagonal link)
      for(int k=1; k<=alpha->mylattice->links[i1][2][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][2][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(nx2, ny2, q2, i2); //get the cartesian coordinates for i2

        //setting hopping to one on all third neighbors
        if( rr && q1!=1 )
          t[n][n][i1][i2] = t[n][n][i2][i1] = -pars->xi[2];
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[2];

        if( pars->e2 ) //unit cell doubling
        {
          if( ny1%2==0 ) {
            if( q1==2 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
          } else
            if( q1==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
/*
          if( ((q1==1)&&(q2==2))||((q1==2)&&(q2==1)) ) //base link of dn triangle
            if( ny1%2==0 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }

          if( ((q1==1)&&(q2==0))||((q1==2)&&(q2==0)) ) //left link of up triangle
            if( ny1%2==1 ) {
              t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
            }
*/
        }

        // boundary conditions
        if( pars->apx && ( abs(nx1 - nx2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
        if( pars->apy && ( abs(ny1 - ny2) > L/2 ) ) {
          t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
        }
      } //end loop third neighbors


    } //end loop flavors
  } //end loop sites
}

void u1real::set_hopping3(double *tt1, double *tt2)
{
  u1hybrid::set_hopping3(tt1, tt2, pars->apx);
}

void u1real::print() {}

void u1real::print_avgs()
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

int u1real::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO liquid (id, sites, lattice, txtdsc, NS, N0, apx, apy, mc_length, nbin, e2, TR, gR, ";
  os << "l3, l2, l1, xi1, xi2, xi3, dd1, dd2, dd3, a1, a2, a3, b1, b2, b3, ";
  os << "P1, dP1, P2, dP2, P3, dP3 ) VALUES (";

  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << (int)pars->apx << ", " << (int)pars->apy << ", " << mc_length << ", " << nk << ", ";
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

