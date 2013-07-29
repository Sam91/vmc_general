#include "u1hybthree.h"
#include "spinone.h"
#include <iomanip>

//define a hybridized unpaired fermionic state on a general lattice for NS flavors.

u1hybthree::u1hybthree(int n, int l) : u1hybrid( n, l )
{
  pars = new parameters;
  bpars = pars;


  pars->N = N;

  pars->ap = true;
  
  pars->t1 = 0.;
  pars->t2 = 0.;
  pars->r2 = 0.;
  pars->phi1=0.;pars->eta=0.;pars->phi=0.;pars->theta=0.;pars->psi=0.;
  pars->td1=0.;pars->td2=0.;
}

u1hybthree::~u1hybthree()
{
  //delete[] pars->NF;
  //delete[] pars->mu;
  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1hybthree::set_hopping(double tt1, double tt2)
{
  pars->t1 = tt1; pars->t2 = tt2;
  u1hybrid::set_hopping(tt1, tt2, pars->ap);
}

//set additional flavor-diagonal hoppings on the diagonals of the square lattice (breaking rotation)
void u1hybthree::set_td(double td1, double td2)
{
  cout << "u1hybthree::set_td(" << td1 << "," << td2 << ")\n";
  pars->td1 = td1; pars->td2 = td2;

  int nx2, ny2, i2a, i2b;

  int i1=0;  

  for(int nx=0; nx<L; nx++)
  {
    for(int ny=0; ny<L; ny++)
    {
      //x+y
      nx2 = nx+1; ny2 = ny+1; alpha->mylattice->torus(nx2, ny2); i2a = alpha->mylattice->j(nx2,ny2,0);

      //x-y
      nx2 = nx-1; ny2 = ny+1; alpha->mylattice->torus(nx2, ny2); i2b = alpha->mylattice->j(nx2,ny2,0);

      for(int n=0; n<NS; n++)
      {
        t[n][n][i1][i2a] += td1;
        t[n][n][i1][i2b] += td2;
      }
      i1++;
    }
  }
}

void u1hybthree::set_r2(double r2) //set an additional hopping on the right link of the nn checkerboard
{
  cout << "u1hybthree::set_r2(" << r2 << ")\n";
  pars->r2 = r2;

  int nx, ny, nx2, ny2, i2a;

  for(int i1=0; i1<N; i1++)
  {
    nx = i1 % L;
    ny = i1 / L;

    if( (nx+ny)%2 != 0) continue;

    //x+y
    nx2 = nx+1; ny2 = ny+1; alpha->mylattice->torus(nx2, ny2); i2a = alpha->mylattice->j(nx2,ny2,0);

    for(int n=0; n<NS; n++)
      t[n][n][i1][i2a] += r2;
  }
}


//set a four-sublattice umbrella state
void u1hybthree::set_umb(double eta)
{
  if( (L%4) != 0 ) {
    cout << "ERROR: here, we want a 4-sublattice order\n";
    exit(-1);
  }

  cout << "u1hybthree::set_umb(" << eta << ")\n";
  pars->eta = eta;
  pars->desc = "umbrella";

  int i1=0;
  double ce = cos(eta*M_PI);
  double se = sin(eta*M_PI);

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {
      d[0][i1] = ce;

      if( (i%2)==0 && (j%2)==0 ) {
        d[1][i1] = se; d[2][i1] = 0.;
      }
      if( (i%2)==1 && (j%2)==1 ) {
        d[1][i1] = 0.; d[2][i1] = -se;
      }
      if( (i%2)==1 && (j%2)==0 ) {
        d[1][i1] = 0.; d[2][i1] = se;
      }
      if( (i%2)==0 && (j%2)==1 ) {
        d[1][i1] = -se; d[2][i1] = 0.;
      }

      i1++;
    }
  }
}

//set a three-sublattice umbrella state
void u1hybthree::set_three(double eta)
{
  if( (L%3) != 0 ) {
    cout << "ERROR: here, we want a 3-sublattice order\n";
    exit(-1);
  }

  cout << "u1hybthree::set_three(" << eta << ")\n";
  pars->eta = eta;
  pars->desc = "three";

  int i1=0;
  double ce = cos(eta*M_PI); //perpendicular directors when cos^2 = 1/3
  double se = sin(eta*M_PI);

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {
      d[0][i1] = ce;

      switch( ((i%3)+(j%3))%3 )
      {
        case 0: //A-sites
          d[1][i1] = se;
          d[2][i1] = 0.;
          break;

        case 1: //B-sites
          d[1][i1] = -se/2.;
          d[2][i1] = -sqrt(3.)*se/2.;
          break;

        case 2: //C-sites
          d[1][i1] = -se/2.;
          d[2][i1] =  sqrt(3.)*se/2.;
          break;
      }

      i1++;
    }
  }
}

//set a general (real) three-sublattice state
void u1hybthree::set_three(double phi, double theta, double psi)
{
  if( (L%3) != 0 ) {
    cout << "ERROR: here, we want a 3-sublattice order\n";
    exit(-1);
  }

  cout << "u1hybthree::set_three(" << phi << ", " << theta << ", " << psi << ")\n";
  pars->phi = phi;
  pars->theta = theta;
  pars->psi = psi;
  pars->desc = "three-general";

  int i1=0;
  double ce1 = cos(phi*M_PI); //perpendicular directors when cos^2 = 1/3
  double se1 = sin(phi*M_PI);
  double ce2 = cos(theta*M_PI);
  double se2 = sin(theta*M_PI);
  double ce3 = cos(psi*M_PI);
  double se3 = sin(psi*M_PI);

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {

      switch( ((i%3)+(j%3))%3 )
      {
        case 0: //A-sites
          d[0][i1] = 1.;
          d[1][i1] = 0.;
          d[2][i1] = 0.;
          break;

        case 1: //B-sites
          d[0][i1] = se1;
          d[1][i1] = ce1;
          d[2][i1] = 0.;
          break;

        case 2: //C-sites
          d[0][i1] = ce2;
          d[1][i1] = se2*se3;
          d[2][i1] = se2*ce3;
          break;
      }

      i1++;
    }
  }
}

//set a general (complex) three sublattice state
void u1hybthree::set_three(double phi1, double eta, double phi, double theta, double psi)
{
#if !WFC
  cout << "ERROR: compile with WFC=true\n";
  exit(-1);
#else

  if( (L%3) != 0 ) {
    cout << "ERROR: here, we want a 3-sublattice order\n";
    exit(-1);
  }

  cout << "u1hybthree::set_three(phi1=" << phi1 << ",eta=" << eta << ",phi=" << phi << ",theta=" << theta << ",psi=" << psi << ")\n";
  pars->phi1 = phi1; //angle on B site (real)
  pars->eta = eta;
  pars->phi = phi;
  pars->theta = theta;
  pars->psi = psi;
  pars->desc = "three-general";

  spinone* sp1 = new spinone(eta*M_PI, phi*M_PI, theta*M_PI, psi*M_PI);

  int i1=0;
  double ce1 = cos(phi1*M_PI); //perpendicular directors when cos^2 = 1/3
  double se1 = sin(phi1*M_PI);
  double ce2 = cos(theta*M_PI);
  double se2 = sin(theta*M_PI);
  double ce3 = cos(psi*M_PI);
  double se3 = sin(psi*M_PI);

  cout << "C: " << complex<double>(sp1->getU(2), sp1->getV(2)) << "; " << complex<double>(sp1->getU(0), sp1->getV(0)) << "; " << complex<double>(sp1->getU(1), sp1->getV(1)) << "\n";
  cout << "r: " << ce2 << "; " << se2*se3 << "; " << se2*ce3 << "\n";

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {

      switch( ((i%3)+(j%3))%3 )
      {
        case 0: //A-sites
          d[0][i1] = 1.;
          d[1][i1] = 0.;
          d[2][i1] = 0.;
          break;

        case 1: //B-sites
          d[0][i1] = se1;
          d[1][i1] = ce1;
          d[2][i1] = 0.;
          break;

        case 2: //C-sites
          //d[0][i1] = ce2;
          //d[1][i1] = se2*se3;
          //d[2][i1] = se2*ce3;
          d[0][i1] = complex<double>(sp1->getU(2), sp1->getV(2));
          d[1][i1] = complex<double>(sp1->getU(0), sp1->getV(0));
          d[2][i1] = complex<double>(sp1->getU(1), sp1->getV(1));
          break;
      }

      i1++;
    }
  }
#endif
}

//set a four-sublattice distorted umbrella state
void u1hybthree::set_four(double eta)
{
  if( (L%4) != 0 ) {
    cout << "ERROR: here, we want a 4-sublattice order\n";
    exit(-1);
  }

  cout << "u1hybthree::set_four(" << eta << ")\n";
  pars->eta = eta;
  pars->desc = "four";

  //alpha->mylattice->set_triangular();
  //alpha->mylattice->set_checkerboard();
  //alpha->mylattice->print();

  int i1=0;
  double ce = cos(eta*M_PI);
  double se = sin(eta*M_PI);

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++)
    {
      d[0][i1] = ce;

      if( (i%2)==0 && (j%2)==0 ) {
        d[1][i1] =  se; d[2][i1] = 0.;
      }
      if( (i%2)==1 && (j%2)==1 ) {
        d[1][i1] = -se; d[2][i1] = 0.;
      }
      if( (i%2)==0 && (j%2)==1 ) {
        d[1][i1] = 0.; d[2][i1] =  se;
      }
      if( (i%2)==1 && (j%2)==0 ) {
        d[1][i1] = 0.; d[2][i1] = -se;
      }

      i1++;
    }
  }
}

//set a four-sublattice xxyz state (ik is the state that is twice occupied in a in a 4-site sublattice)
void u1hybthree::set_xxyz(int ik)
{
  if( (L%4) != 0 ) {
    cout << "ERROR: here, we want a 4-sublattice order\n";
    exit(-1);
  }

  cout << "huseelser::set_xxyz(" << ik << ")\n";
  pars->eta = 0;
  stringstream s; s << "xxyz-r-" << ik;
  pars->desc = s.str();

  int i1=0;

  //determine the other two states
  int k1, k2; k1=0; k2=0;
  for(int k=0; k<NS; k++)
  {
    if( k!= ik ) {
      k1 = k; break;
    }
  }
  for(int k=0; k<NS; k++)
  {
    if( k!= ik && k!= k1) {
      k2 = k; break;
    }
  }

  for(int i=0; i<L; i++)
  {
    for(int j=0; j<L; j++)
    {
      //d[ik][i1] = 1./sqrt(3.);

      if( (i%2)==0 && (j%2)==0 ) {
        d[ik][i1] = 1.;
        d[k1][i1] = 0.;
        //d[k2][i1] = sqrt(2./3.);
        d[k2][i1] = 0.;
      }
      if( (i%2)==1 && (j%2)==1 ) {
        d[ik][i1] = 1.;
        d[k1][i1] = 0.;
        //d[k2][i1] = sqrt(2./3.);
        d[k2][i1] = 0.;
      }
      if( (i%2)==0 && (j%2)==1 ) {
        d[ik][i1] = 0.;
        d[k1][i1] =  sqrt(.5);
        //d[k2][i1] = -sqrt(1./6.);
        d[k2][i1] = sqrt(.5);
      }
      if( (i%2)==1 && (j%2)==0 ) {
        d[ik][i1] = 0.;
        d[k1][i1] = -sqrt(.5);
        //d[k2][i1] = -sqrt(1./6.);
        d[k2][i1] = sqrt(.5);
      }

      i1++;
    }
  }
}

void u1hybthree::print()
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

/*void u1hybthree::print_t(int n1, int n2)
{
  int i1=0;
  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++) {
      cout << t[n1][n2][i1][] << ", ";
      i1++;
    }
    cout << "\n";
  }
}*/

int u1hybthree::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO u1hybthree (id, sites, lattice, NS, N0, t1, t2, h, phi1, eta, phi, theta, psi, muz, ap, txtdsc, td1, td2, r2, ";
  os << "P1, dP1, P2, dP2) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', ";
  os << NS << ", " << NF[0] << ", ";
  os << "round(" << pars->t1 << ",3), round(" << pars->t2 << ",3), round(" << h << ",3),";
  os << "round(" << pars->phi1 << ",3), round(" << pars->eta << ",3), round(" << pars->phi << ", 3), round(" << pars->theta << ",3), round(" << pars->psi << ", 3), " << mu[0] << ", ";
  os << pars->ap << ", '" << pars->desc << "', round(" << pars->td1 << ",3), round(" << pars->td2 << ",3), round(" << pars->r2 << ",3), ";
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

