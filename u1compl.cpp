#include "u1compl.h"
#include <iomanip>

// Define an unpaired real fermionic state on the Kagome lattice with first, second, and diagonal neighbor hopping
// The state may have unit-vell doubling (epsilon) and an additional rotation staggering of the hopping real (tau_R)

u1compl::u1compl(int l, int q) : u1real( l, q )
{}

u1compl::~u1compl()
{}

//set a complexe t1-t2-t3 hopping state on Kagome 
void u1compl::set_hoppingk(double mu0)
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

  double srr; if( rr ) srr = -1.; else srr = 1.;
  double sri; if( pars->TR ) sri = -1.; else sri = 1.;

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
        if( (q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1) )
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[0]*complex<double>( srr*cos(M_PI*pars->a[0]), sri*sin(M_PI*pars->a[0]) );
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[0]*complex<double>( cos(M_PI*pars->a[0]), sin(M_PI*pars->a[0]) );

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
        if( (q1==0 && k==1)||(q1==1 && k==2)||(q1==2 && k==1) )
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[1]*complex<double>( srr*cos(M_PI*pars->a[1]), sri*sin(M_PI*pars->a[1]) );
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[1]*complex<double>( cos(M_PI*pars->a[1]), sin(M_PI*pars->a[1]) );

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
        if( q1!=1 )
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[2]*complex<double>( srr*cos(M_PI*pars->a[2]), sri*sin(M_PI*pars->a[2]) );
        else
          t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[2]*complex<double>( cos(M_PI*pars->a[2]), sin(M_PI*pars->a[2]) );

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

