#include "u1dirac.h"
#include "spinone.h"
#include <iomanip>

// Define an unpaired real fermionic state on the Kagome lattice with first, second, and diagonal neighbor hopping
// The state may have unit-cell doubling (epsilon) and an additional rotation staggering of the hopping; (taur, taui)

u1dirac::u1dirac(int l, int q) : u1hybrid( l, q )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N;

  pars->apx = false;
  pars->apy = false;
}

u1dirac::~u1dirac()
{
  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1dirac::set_hopping()
{
  u1hybrid::set_hopping(1., 0, pars->apx);
}

void u1dirac::set_hoppingk()
{
  for(int n=0; n<NS; n++) mu[n] = 0.;

  int i2;
  int nx2, ny2, q2, nx1, ny1, q1;

  //set the hopping matrix to zero
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n1=0; n1<NS; n1++)
        for(int n2=0; n2<NS; n2++) t[n1][n2][i][j] = 0.;

  for(int i1=0; i1<N; i1++) //loop over all sites
  {
    alpha->mylattice->getnq(nx1, ny1, q1, i1); //get the cartesian coordinates for i1.

    for(int n=0; n<NS; n++) //loop over all flavors (2)
    {
      for(int k=1; k<=alpha->mylattice->links[i1][0][0]; k++) //loop over nearest neighbor sites
      {
        i2 = alpha->mylattice->links[i1][0][k];

        t[n][n][i1][i2] = t[n][n][i2][i1] = -1.; //setting hopping to one on all nearest neighbors

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(nx2, ny2, q2, i2); //get the cartesian coordinates for i2

        if( (q1==0) && (q2==1) ) //base link of up triangle
          if( nx1%2==1 )
            t[n][n][i1][i2] = t[n][n][i2][i1] = 1.;

        if( (q1==0) && (q2==2) ) //left link of up triangle
          if( nx1%2==0 )
            t[n][n][i1][i2] = t[n][n][i2][i1] = 1.;

        if( (q1==1) && (q2==2) ) //right link of up triangle
            t[n][n][i1][i2] = t[n][n][i2][i1] = 1.;


        // boundary conditions
        if( pars->apx && ( abs(nx1 - nx2) > L/2 ) ) {t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;}
        if( pars->apy && ( abs(ny1 - ny2) > L/2 ) ) {t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;}
      }

    }  
  } //end loop over all sites
}

void u1dirac::set_hopping3(double *tt1, double *tt2)
{
  u1hybrid::set_hopping3(tt1, tt2, pars->ap);
}

void u1dirac::print() {}

void u1dirac::print_avgs()
{
  cout<< "ff={";
  for(int no=0; no<NO-1; no++) {
    cout << std::fixed << setprecision(4) << average[no] << ", ";
  }
  cout << std::fixed << setprecision(4) << average[NO-1] << "};"<<endl;
  cout<< "ffsigma={";
  for(int no=0; no<NO-1; no++) {
    cout << std::fixed << setprecision(4) << sigma[no] << ", ";
  }
  cout << std::fixed << setprecision(4) << sigma[NO-1] << "};"<<endl;
}


int u1dirac::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, ap, t1, t2, h, t1b, t2b, t1c, t2c, phi1, phi2, mu, j1a, j1b, j1c, mc_length, nk, ";
  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4, P1p, dP1p) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << (int)pars->ap << ", ";
  os << "round(" << js[0] << ",3), round(" << js[3] << ",3), round(" << js[5] << ",3), ";
  os << mc_length << ", " << nk << ", ";
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ", " << average[4] << ", " << sigma[4];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

