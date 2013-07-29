#include "matrix.h"
#include "helperfunctions.h"

int main()
{
//  test_diag(3500, 5);


  complex<double> **m = createcomplex(3);
  complex<double> **m1 = createcomplex(3);

  m[0][0] = .1;
  m[0][1] = -.1;
  m[0][2] = .5;
  m[1][0] = 5.6;
  m[1][1] = -.3;
  m[1][2] = .2;
  m[2][0] = 8.7;
  m[2][1] = .7;
  m[2][2] = -.5;

  for(int i=0; i<3; i++)
    for(int j=i; j<3; j++)
      m[i][j] = (m[i][j]+m[j][i])/2.;
  for(int i=0; i<3; i++)
    for(int j=0; j<i; j++)
      m[i][j] = m[j][i];

  write_m(m, 3);

  double *v = new double[3];

  complex<double> *inv = new complex<double>[9];
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) inv[i*3+j] = m[i][j];

  complex<double> det1 = inverse(m, m1, 3);
  complex<double> det2 = inverse(inv, 3);

  cout << "det1: " << det1 << "\n";
  write_m(m1, 3);

  //det = inverse_m(m, m1, 3);

  cout << "det2: " << det2 << "\n";
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) m1[i][j] = inv[i*3+j];
  write_m(m1, 3);
/*
  eigvects(m, v, 3);

  cout << "eigv:\n";
  write_m(m, 3);

  stringstream s; s << "hello ";
  cout << s.str() << "\n";
  s << "world " << 3;
  cout << s.str() << "\n";
*/
  return 0;
}
