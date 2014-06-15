#include "u1two.h"
#include "matrix.h"

//define a determinental U(1) state for TWO flavors with equal filling (wf is a single determinant).

u1two::u1two(int l) : twoflavor( l )
{
  if(NS!=2) {
    cout << "Must be compiled with NS=2\n";
    exit(-2);
  }
  
  if(N%NS != 0) {
    cout << "N (number of sites) must be a multiple of " << NS << "\n";
    exit(-1);
  }

  pars = new parameters;
  bpars = pars;
  pars->t = new (double(**[NS]));

  for(int n=0; n<NS; n++)
  {
    pars->t[n] = createdouble( N );
  }

  pars->t1 = 1.;
  pars->t2 = 1.;
  pars->ap = true;
  pars->N = N;
}

u1two::~u1two()
{
  for(int n=0; n<NS; n++)
  {
    destroy(pars->t[n], N);
  }

  delete[] pars->t;
  bpars = NULL;
  delete pars;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1two::set_hopping(double t, double tt)
{
  pars->t1 = t;
  pars->t2 = tt;

  int i2;

  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n=0; n<NS; n++)
        pars->t[n][i][j] = 0.;

  // the labelling of sites is specified in mylattices; it is ix*L+iy.

  for(int i1=0; i1<N; i1++) {
    for(int n=0; n<NS; n++) {
      for(int k=1; k<=alpha->mylattice->connectivity[i1][0][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][0][k];
        if( pars->ap && ( abs(i1%L - i2%L)>L/2 ) ) pars->t[n][i1][i2] += t; //FIXME: not sure if BC work for lattices with more than one site per unit cell
        else pars->t[n][i1][i2] += -t;
      }
      for(int k=1; k<=alpha->mylattice->connectivity[i1][1][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][1][k];
        if( pars->ap && ( abs(i1%L - i2%L)>L/2 ) ) pars->t[n][i1][i2] -= tt; //FIXME: same as above
        else pars->t[n][i1][i2] += tt;
      }
    }
  }
}

void u1two::print()
{
  cout << "wf pars->\n";
  cout << " J1, J2 = " << js[0] << ", " << js[1] << "\n";

  //cout << "t: " << t << "\n";
}

//initialize the U(1) state
void u1two::create_ad()
{
  pars->desc = "U(1)2";

  double **m = createdouble( N );
  double *v = new double[ N ];
  
  double ***adxn = new (double(**[2]));
  adxn[0] = createdouble( N );
  adxn[1] = createdouble( N );

  for(int n=0; n<NS; n++)
  {
    copy_m(m, pars->t[n], N);
    eigvects(m, v, N);
    cout << "v[" << n << "]: ";
    for(int i=0; i<N; i++) cout << v[i] << "; ";
    cout << "\n";

    //take the lowest N2 eigenvectors
    int i;
    for(i=0; i<N2; i++)
      for(int j=0; j<N; j++) adxn[n][j][i] = m[i][j]; //note that the eigenvalue index is the last index in adx!

    if(i<N) {
      //cout << "abs diff " << abs(v[i]-v[i-1]) << "\n";
      if( abs(v[i]-v[i-1]) < 1e-8 )
        cout << "WARN: degenerate state\n";
      cout << "Fermi energy is " << v[i-1] << ", " << v[i] << " at position " << i-1 << "\n"; 
    }
  }
  
  for(int j=0; j<N*N; j++) adx[j] = 0.;
  
  for(int i=0; i<N; i++) //mutliply the two matrices to get adx
    for(int j=0; j<N; j++)
      for(int n=0; n<N2; n++)
        adx[i*N+j] += adxn[0][i][n]*adxn[1][j][n];

  destroy(m, N);
  delete[] v;
  destroy(adxn[0], N);
  destroy(adxn[1], N);
  delete[] adxn;

  normalize( 4. );
}
int u1two::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, ap, ";
  os << "t1, t2, P1, dP1, P2, dP2) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', " << NS << ", " << N2 << ", " << (int)pars->ap << ", '" << pars->desc << "'";
  os << ", round(" << pars->t1 << ",3), round(" << pars->t2 << ",3), ";
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

