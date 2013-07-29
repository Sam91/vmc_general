#include "u1hybtwo.h"
#include "spinone.h"
#include <iomanip>

//define a hybridized unpaired fermionic state on a general lattice for NS flavors.

u1hybtwo::u1hybtwo(int n, int l) : u1hybrid( n, l )
{
  pars = new parameters;
  bpars = pars;

  pars->N = N;

  pars->ap = true;
  
  pars->t1  = 0.;
  pars->t1b = 0.;
  pars->t1c = 0.;
  pars->t2  = 0.;
  pars->t2b = 0.;
  pars->t2c = 0.;

  pars->phi1=0.;pars->phi2=0.;
}

u1hybtwo::~u1hybtwo()
{
  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1hybtwo::set_hopping(double tt1, double tt2)
{
  pars->t1 = tt1; pars->t2 = tt2;
  pars->t1b = tt1; pars->t2b = tt2;
  pars->t1c = tt1; pars->t2c = tt2;
  u1hybrid::set_hopping(tt1, tt2, pars->ap);
}

void u1hybtwo::set_hopping3(double *tt1, double *tt2)
{
  pars->t1 = tt1[0];
  pars->t1b = tt1[1];
  pars->t1c = tt1[2];
  pars->t2 = tt2[0];
  pars->t2b = tt2[1];
  pars->t2c = tt2[2];
  u1hybrid::set_hopping3(tt1, tt2, pars->ap);
}

//set a spiral state
void u1hybtwo::set_spiral(double phi1, double phi2)
{
  cout << "u1hybtwo::set_spiral(" << phi1 << ", " << phi2 << ")\n";
  pars->phi1 = phi1;
  pars->phi2 = phi2;
  pars->desc = "spiral";

  int i1=0;
  double cp, sp;

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {
      cp = cos( (i*phi1 + j*phi2)*tPi/2. );
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

void u1hybtwo::print()
{
  cout << "wf pars->\n";
  cout << " J1, J2 = " << js[0] << ", " << js[1] << "\n";

  cout << " d[n][i] =\n";
  for(int i=0; i<N; i++) {
    cout << " " << std::fixed << setprecision(3) << i << " [" << d[0][i] << ", " << d[1][i] << ", " << d[2][i] << "]\n";
  }
  //cout << "t: " << t << "\n";
  cout << "NF: ";
  for(int n=0; n<NS; n++) cout << NF[n] << ", ";
  cout << "\n";
}

int u1hybtwo::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, ap, t1, t2, h, t1b, t2b, t1c, t2c, phi1, phi2, mu, j1a, j1b, j1c, mc_length, nk, ";
  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4, P1p, dP1p) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << NF[0] << ", " << (int)pars->ap << ", ";
  os << "round(" << pars->t1 << ",3), round(" << pars->t2 << ",3), round(" << h << ",3), ";
  os << "round(" << pars->t1b << ",3), round(" << pars->t2b << ",3), ";
  os << "round(" << pars->t1c << ",3), round(" << pars->t2c << ",3), ";
  os << "round(" << pars->phi1 << ",3), round(" << pars->phi2 << ",3), round(" << mu[0] << "), ";
  os << "round(" << js[0] << ",3), round(" << js[3] << ",3), round(" << js[5] << ",3), ";
  os << mc_length << ", " << nk << ", ";
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ", " << average[4] << ", " << sigma[4];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

