#include "u1real.h"
#include <iomanip>

// Define an unpaired real fermionic state on the Kagome lattice with first, second, and diagonal neighbor hopping
// The state may have unit-vell doubling (epsilon) and an additional rotation staggering of the hopping real (tau_R)

u1real::u1real( int l ) : u1hybrid( l )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N;

  pars->ap = new bool[DIM];

  for(int d=0; d<DIM; d++) pars->ap[d] = false;

  int nxi = 3;
  pars->xi = new double[nxi]; for(int i=0; i<nxi; i++) pars->xi[i] = 0.;
}

u1real::~u1real()
{
  delete[] pars->xi;
  delete[] pars->ap;
  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1real::set_hopping()
{
  u1hybrid::set_hopping(1., 0, pars->ap[0]);
}

//set the real hopping matrix on three neighbors
void u1real::set_hoppingk(double mu0)
{
  for(int n=0; n<NS; n++) mu[n] = mu0;

  int i2;
  int q1, q2;

  cout << "ap = [";
  for(int d=0; d<DIM; d++) cout << pars->ap[d] << "; ";
  cout  << "]" << endl;

  cout << "xi = ["<< pars->xi[0] << "; " << pars->xi[1] << "; " << pars->xi[2] << "]" << endl;

  //set the hopping matrix to zero
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int nn1=0; nn1<NS; nn1++)
        for(int nn2=0; nn2<NS; nn2++) t[nn1][nn2][i][j] = 0.;

  int *n1 = new int[DIM];
  int *n2 = new int[DIM];

  for(int i1=0; i1<N; i1++) //loop over all sites
  {
    alpha->mylattice->getnq(n1, q1, i1); //get the cartesian coordinates for i1.

    for(int n=0; n<NS; n++) //loop over flavors (two)
    {
      //loop over first neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][0][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][0][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); // get the cartesian coordinates for i2

        // set hopping on nearest neighbors
        t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[0];

        // boundary conditions
        for(int d=0; d<DIM; d++)
        {
          if( pars->ap[d] && ( abs(n1[d] - n2[d]) > L/2 ) ) {
            t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }
        }
      } //end loop first neighbors

      //loop over second neighbors
      for(int k=1; k<=alpha->mylattice->links[i1][1][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][1][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); //get the cartesian coordinates for i2

        //set hopping on second neighbor according
        t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[1];

        // boundary conditions
        for(int d=0; d<DIM; d++)
        {
          if( pars->ap[d] && ( abs(n1[d] - n2[d]) > L/2 ) ) {
            t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }
        }
      } //end loop second neighbors

      //loop over third neighbors (diagonal link)
      for(int k=1; k<=alpha->mylattice->links[i1][2][0]; k++)
      {
        i2 = alpha->mylattice->links[i1][2][k];

        //cout<<"i1="<<i1<<"; i2="<<i2<<"; k="<<k<<endl;

        alpha->mylattice->getnq(n2, q2, i2); //get the cartesian coordinates for i2

        //setting hopping to one on all third neighbors
        t[n][n][i1][i2] = t[n][n][i2][i1] = pars->xi[2];

        // boundary conditions
        for(int d=0; d<DIM; d++)
        {
          if( pars->ap[d] && ( abs(n1[d] - n2[d]) > L/2 ) ) {
            t[n][n][i1][i2] *= -1.; t[n][n][i2][i1] *= -1.;
          }
        }
      } //end loop third neighbors

    } //end loop flavors
  } //end loop sites

  delete[] n1; delete[] n2;
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

  os << "INSERT INTO liquid (id, sites, lattice, txtdsc, NS, N0, apx, mc_length, nbin, ";
  os << "l3, xi1, xi2, xi3, ";
  os << "P1, dP1, P2, dP2, P3, dP3 ) VALUES (";

  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << (int)pars->ap[0] << ", " << mc_length << ", " << nk << ", ";

  os << "round(" << mu[0] << ",3), ";
  os << "round(" << pars->xi[0] << ",3), round(" << pars->xi[1] << ",3), round(" << pars->xi[2] << ",3), ";

  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

//int u1real::insert_file(const char *name)
int u1real::insert_file(string name)
{
  int kk = 4+2*NO;
  double* ftmp = new double[ kk ];
  //cout << "allocating: " << kk << "\n";

  int idx = 0;
  
  ftmp[idx] = N; idx++;
  for(int i=0; i<3 ; i++) { ftmp[i+idx] = pars->xi[i];} idx += 3;
  for(int i=0; i<NO; i++) { ftmp[i+idx] = average[i]; } idx += NO;
  for(int i=0; i<NO; i++) { ftmp[i+idx] = sigma[i];   } idx += NO;

  int r = fappend(ftmp, kk, name);

  delete[] ftmp;
  ftmp = nullptr;

  return r;
}

