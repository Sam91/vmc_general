/*
 * Some functions specific to the implementation of ordered states on the triangular lattice
 *
**/

#include "spinone.h"
#include "matrix.h"

void spinone::init()
{
  //cout << "Creating the spinone data\n";
  u = new double[3];
  v = new double[3];

  D = createcomplex(3);
}

spinone::spinone()
{
  init();

  u[0] = 0.; u[1] = 0.; u[2] = 1.;
  v[0] = 0.; v[1] = 0.; v[2] = 0.;
  createD();
}

spinone::~spinone()
{
  //cout << "Destroying spinone\n";
  delete[] u; delete[] v;
  destroy(D, 3);
}

//Spin-one state from u and v vectors
spinone::spinone(double *u0, double *v0)
{
  init();

  double n = 0.;
  double m = 0.;

  for(int i=0; i<3; i++){
    u[i] = u0[i];
    v[i] = v0[i];
    n += u0[i]*u0[i]+v0[i]*v0[i];
    m += u0[i]*v0[i];
  }

  if( abs(n-1.)>TINY || abs(m)>TINY )
    cout << "WARN: Those u and v do not correspond to physical spin states!\n";

  //create the D-matrix
  createD();
}

//create a spin-one state from eta (length of u) and Euler angles
spinone::spinone(double eta, double phi, double theta, double psi)
{
  init();

  //starting values of u and v vectors
  u[0] = 0.; u[1] = 0.; u[2] = cos(eta);
  v[0] = 0.; v[1] = sin(eta); v[2] = 0;

  //rotate by euler angles
  rotate(phi, theta, psi);

  //create the D-matrix from the rotated u and v vectors
  createD();
}

//copy constructor
spinone::spinone(const spinone &s0)
{
  init();
  cout << "Executing copy constructor\n";
  for(int i=0; i<3; i++) {
    u[i] = s0.u[i];
    v[i] = s0.v[i];
    for(int j=0; j<3; j++) D[i][j] = s0.D[i][j];
  }
}

void spinone::operator = (const spinone &s0)
{
  cout << "spinone: Executing assignment operation\n";
  for(int i=0; i<3; i++) {
    u[i] = s0.u[i];
    v[i] = s0.v[i];
    for(int j=0; j<3; j++) D[i][j] = s0.D[i][j];
  }
}

double spinone::getU(int i)
{
  if(i<0||i>2) {
    cout << "getU invalid parameter\n";
    return 0.;
  }
  return u[i];
}
double spinone::getV(int i)
{
  if(i<0||i>2) {
    cout << "getV invalid parameter\n";
    return 0.;
  }
  return v[i];
}

//Rotate the u and v vectors by the Euler angles
void spinone::rotate(double phi, double theta, double psi)
{
  //first, create the rotation matrix (3-1-3)
  double **r = createdouble(3);
  double **m1 = createdouble(3);

  //phi rotation around z
  r[0][0] =  cos(phi); r[0][1] = sin(phi); r[0][2] = 0.;
  r[1][0] = -sin(phi); r[1][1] = cos(phi); r[1][2] = 0.;
  r[2][0] = 0.; r[2][1] = 0.; r[2][2] = 1.;

  //theta rotation around x
  m1[0][0] = 1.; m1[0][1] = 0.; m1[0][2] = 0.;
  m1[1][0] = 0.; m1[1][1] =  cos(theta); m1[1][2] = sin(theta);
  m1[2][0] = 0.; m1[2][1] = -sin(theta); m1[2][2] = cos(theta);

  mult(r, m1, 3);

  //phi rotation around z
  m1[0][0] =  cos(psi); m1[0][1] = sin(psi); m1[0][2] = 0.;
  m1[1][0] = -sin(psi); m1[1][1] = cos(psi); m1[1][2] = 0.;
  m1[2][0] = 0.; m1[2][1] = 0.; m1[2][2] = 1.;

  mult(r, m1, 3);

  //rotate the vectors u and v with the r-matrix
  mult(u, r, 3);
  mult(v, r, 3);

  destroy(r,3); destroy(m1,3);
}

//print out the data for this state
void spinone::print()
{
  cout << "U vector: (";
  for(int i=0; i<3; i++) cout << u[i] << ", ";
  cout << ")\n";

  cout << "V vector: (";
  for(int i=0; i<3; i++) cout << v[i] << ", ";
  cout << ")\n";

  cout << "D-matrix:\n";
  write_m(D, 3);
}

//create the D matrix from the u and v vectors
void spinone::createD()
{
  //write down the D matrix
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      D[i][j] = complex<double>(u[i]*u[j]+v[i]*v[j], u[i]*v[j]-v[i]*u[j]);

  //print();
}

bool spinone::uncoupled_z(){ 
  return (abs(u[2])<TINY && abs(v[2])<TINY );
}

// calculate the projected ground state or highest excited state of the fermionic operator F = D_ab f_a^dagger f_b
// F is block diagonal: the particle block ist simply D, the hole block is 1-D*
// Therefore, we only need to diagonalize D
void spinone::projectD()
{
  complex<double> **m2 = createcomplex(3);
  double *en = new double[3];

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++)
      m2[i][j] = D[i][j];
  }

  eigvects(m2, en, 3);
  cout << "Spectrum (particle rep):\n";
  for(int i=0; i<3; i++) {
    normalize(m2[i], 3);
    cout << "en = " << (abs(en[i])<TINY?0:en[i]) << ": [";
    for(int j=0; j<3; j++) cout << (abs(m2[i][j])<TINY?0:m2[i][j]) << (j==2?"]\n":", ");
  }
  cout << "\n";

  delete[] en; destroy(m2,3);
}


/* The following code is not used anymore because I know now what the eigenvalues and eigenvectors of the F matrix are
void spinone::projectD()
{
  //Fermionic operator is an 8x8 matrix
  complex<double> **m2 = createcomplex(8);
  double *e = new double[8];

  //Build up the Fermionic operator in the basis (0,x,y,z,xy,yz,zx,xyz)
  for(int i=0; i<8; i++) {
    for(int j=0; j<8; j++)
      m2[i][j] = 0.;
  }
  for(int i=1; i<4; i++) m2[i][i] = D[i-1][i-1];

  m2[4][4] = D[1][1]+D[2][2];
  m2[5][5] = D[2][2]+D[0][0];
  m2[6][6] = D[0][0]+D[1][1];

  m2[7][7] = D[0][0]+D[1][1]+D[2][2];

  m2[1][2] = D[0][1]; m2[2][1] = D[1][0];
  m2[2][3] = D[1][2]; m2[3][2] = D[2][1];
  m2[3][1] = D[2][0]; m2[1][3] = D[0][2];

  m2[4][5] = -D[1][0]; m2[5][4] = -D[0][1];
  m2[5][6] = -D[2][1]; m2[6][5] = -D[1][2];
  m2[6][4] = -D[0][2]; m2[4][6] = -D[2][0];

  cout << "D-operator in fermion space:\n";
  write_m(m2, 8);

  //diagonalize this matrix
  eigvects(m2, e, 8);

  double e0 = e[0];
  bool printit;
  complex<double> *vect = new complex<double>[3];

  cout << "spectrum: ";
  for(int i=0; i<8; i++) cout << (abs(e[i])<TINY?0:e[i]) << ", ";
  cout << "\n";

  cout << "Lowest projected eigenvectors (e=" << e0 << "):\n";
  for(int i=0; i<8; i++) {
    if( abs(e[i]-e0)<TINY ) {
      printit = false;
      for(int j=1; j<4; j++) {
        if( abs(m2[i][j])>TINY) printit = true;
        vect[j-1] = m2[i][j];
      }
      if(printit) {
        normalize(vect, 3);
        cout << "particle: ";
        for(int j=0; j<3; j++) cout << ((abs(vect[j])<TINY)?0:vect[j]) << ", ";
        cout << "\n";
      }
      printit = false;
      for(int j=4; j<7; j++) {
        if( abs(m2[i][j])>TINY) printit = true;
        vect[j-4] = m2[i][j];
      }
      if(printit) {
        normalize(vect, 3);
        cout << "hole: ";
        for(int j=0; j<3; j++) cout << ((abs(vect[j])<TINY)?0:vect[j]) << ", ";
        cout << "\n";
      }
    }
  }
  cout << "\n";

  e0 = e[7];
  cout << "Highest projected eigenvectors (e=" << e0 << "):\n";
  for(int i=0; i<8; i++) {
    printit = false;
    if( abs(e0-e[i])<TINY ) {
      for(int j=1; j<4; j++) {
        if( abs(m2[i][j])>TINY) printit = true;
        vect[j-1] = m2[i][j];
      }
      if(printit) {
        normalize(vect, 3);
        cout << "particle: ";
        for(int j=0; j<3; j++) cout << ((abs(vect[j])<TINY)?0:vect[j]) << ", ";
        cout << "\n";
      }
    }
    printit = false;
    if( abs(e0-e[i])<TINY ) {
      for(int j=4; j<7; j++) {
        if( abs(m2[i][j])>TINY) printit = true;
        vect[j-4] = m2[i][j];
      }
      if(printit) {
        normalize(vect, 3);
        cout << "hole: ";
        for(int j=0; j<3; j++) cout << ((abs(vect[j])<TINY)?0:vect[j]) << ", ";
        cout << "\n";
      }
    }
  }
  cout << "\n";

  delete[] e; delete[] vect;
  destroy(m2,8);
}
*/

//normalize a complex vector and fix the phase such that im(vect).re(vect)=0
void spinone::normalize(complex<double> *vect, int n)
{
  double norm = 0.; double uv = 0.; double u2 = 0.; double v2 = 0.;
  for(int i=0; i<n; i++) {
    norm += pow(abs(vect[i]),2);
    uv += vect[i].real()*vect[i].imag();
    u2 += pow(vect[i].real(), 2);
    v2 += pow(vect[i].imag(), 2);
  }
  if(abs(uv)>TINY) {
    double ph = atan(2.*uv/(v2-u2))/2.;
    for(int i=0; i<n; i++) vect[i] = vect[i]*complex<double>(cos(ph)/sqrt(norm),sin(ph)/sqrt(norm));
  } else
    for(int i=0; i<n; i++) vect[i] = vect[i]/sqrt(norm);
}
