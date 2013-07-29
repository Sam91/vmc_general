#include "helperfunctions.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <iomanip>

//Allocates and returns a square 2d array of integers
int** createint(int n) {
  int** ret = new (int(*[n]));
  for(int i=0;i<n;i++) ret[i] = new int[n];
  return ret;
}
int** createint(int n, int m) {
  int** ret = new (int(*[n]));
  for(int i=0;i<n;i++) ret[i] = new int[m];
  return ret;
}

//Allocates and returns a square 2d array of real numbers
double** createdouble(int n) {
  double** ret = new (double(*[n]));
  for(int i=0;i<n;i++) ret[i] = new double[n];
  return ret;
}
double** createdouble(int n, int m) {
  double** ret = new (double(*[n]));
  for(int i=0;i<n;i++) ret[i] = new double[m];
  return ret;
}
complex<double>*** createcomplex(int n, int m, int l) {
  complex<double>*** ret = new (complex<double>(**[n]));
  for(int i=0;i<n;i++) {
    ret[i] = new (complex<double>(*[m]));
    for(int j=0;j<m;j++) ret[i][j] = new complex<double>[l];
  }
  return ret;
}

//Allocates and returns a square 2d array of complex numbers
complex<double>** createcomplex(int n) {
  complex<double>** ret = new (complex<double>(*[n]));
  for(int i=0;i<n;i++) ret[i] = new complex<double>[n];
  return ret;
}
complex<double>** createcomplex(int n, int m) {
  complex<double>** ret = new (complex<double>(*[n]));
  for(int i=0;i<n;i++) ret[i] = new complex<double>[m];
  return ret;
}

//Deallocates the 2d array ptr
void destroy( int** ptr, int n) {
  for(int i=0;i<n;i++) delete[] ptr[i];
  delete[] ptr;
}
void destroy( double** ptr, int n) {
  for(int i=0;i<n;i++) delete[] ptr[i];
  delete[] ptr;
}
void destroy(double*** ptr, int n, int m) {
  for(int i=0;i<n;i++) 
    for(int j=0; j<n; j++)
      delete[] ptr[i][j];
  delete[] ptr;
}
void destroy( complex<double>** ptr, int n) {
  for(int i=0;i<n;i++) delete[] ptr[i];
  delete[] ptr;
}

void write_m(double **a, int n)
{
  cout << std::setprecision(9) << std::showpoint << std::fixed;

  cout << "{";
  for(int i=0;i<n;i++) {
    cout << "{";
    for(int j=0;j<n;j++) {
      cout << (abs(a[i][j])>1e-10||a[i][j]!=a[i][j]?a[i][j]:0);
      if (j<n-1) cout << ", ";
    }
    if (i<n-1) cout << "},\n";
    else cout << "}}\n";
  }
}
void write_m(complex<double> **a, int n)
{
  cout << "{";
  for(int i=0;i<n;i++) {
    cout << "{";
    for(int j=0;j<n;j++) {
      cout << (abs(a[i][j].real())>1e-10?a[i][j].real():0) << "+" << (abs(a[i][j].imag())>1e-10?a[i][j].imag():0) << " I";
      if(j<n-1) cout << ", ";
    }
    if (i<n-1) cout << "},\n";
    else cout << "}}\n";
  }
}
void write_m(int **a, int n) {
  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++)
      cout << a[i][j] << " ";
    cout << ";\n";
  }
}
void write_m(double *a, int n)
{
  cout << std::setprecision(9) << std::showpoint << std::fixed;

  cout << "{";
  for(int i=0;i<n;i++) {
    cout << "{";
    for(int j=0;j<n;j++) {
      cout << (abs(a[i*n+j])>1e-10||a[i*n+j]!=a[i*n+j]?a[i*n+j]:0);
      if (j<n-1) cout << ", ";
    }
    if (i<n-1) cout << "},\n";
    else cout << "}}\n";
  }
}
void write_m(complex<double> *a, int n)
{
  cout << std::setprecision(9) << std::showpoint << std::fixed;

  cout << "{";
  for(int i=0;i<n;i++) {
    cout << "{";
    for(int j=0;j<n;j++) {
      cout << (abs(a[i*n+j])>1e-10||a[i*n+j]!=a[i*n+j]?a[i*n+j]:0);
      if (j<n-1) cout << ", ";
    }
    if (i<n-1) cout << "},\n";
    else cout << "}}\n";
  }
}

//save/load a square matrix in binary format to a file
void save_m(double** a, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "wb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return;}
  for(int i=0; i<n; i++)
    fwrite(a[i], sizeof(double), n, pfile);
  fclose(pfile);
}
void load_m(double** a, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "rb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return;}
  for(int i=0; i<n; i++)
    if(fread(a[i], sizeof(double), n, pfile)){};
  fclose(pfile);
}

//save and load a row vector of doubles
int save_v(double* v, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "wb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return -1;}
  fwrite(v, sizeof(double), n, pfile);
  fclose(pfile);
  return 0;
}
int load_v(double* v, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "rb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return -1;}
  if(fread(v, sizeof(double), n, pfile)){};
  fclose(pfile);
  return 0;
}

//save and load a vector of integers
int save_v(int* v, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "wb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return -1;}
  fwrite(v, sizeof(int), n, pfile);
  fclose(pfile);
  return 0;
}
int load_v(int* v, int n, const char* filename)
{
  FILE* pfile = fopen(filename, "rb");
  if(pfile==NULL) {cout << "Error: cannot open file\n"; return -1;}
  if(fread(v, sizeof(int), n, pfile)){};
  fclose(pfile);
  return 0;
}

/* simple function to write out the remaining time for the calculation */
void estimate_time( double time_start, int to_do, int done )
{
  double delta_time = time(0)-time_start;
  double total_length = (double)to_do/(double)done*delta_time;

  cout << "\nUntil now, we have used " << std::fixed << setprecision(1) << round(delta_time/6.)/10. << " minutes for the runs (" << done << "/" << to_do << " runs)\n";
  cout << "It will take us " << setprecision(2) << round(total_length/36.)/100. << " hours (" << setprecision(1) << round(total_length/6.)/10. << " mins)\n";
  double remaining = (total_length - delta_time)/3600.;
  cout << "Estimate for remaining time: " << setprecision(2) << round(remaining*100.)/100. << " hours (" << setprecision(1) << round(remaining*600.)/10. << " mins)" << endl;
}

bool checkskew(double** a, int n)
{
  //check if the matrix is skew symmetric
  double diff; int counter = 0;
  for(int i=0; i<n; i++) {
    for(int j=0; j<i; j++) {
      if( abs(a[i][j])>1e-6 && abs(a[j][i])>1e-6 ) {
        diff = abs((a[i][j] + a[j][i])/(a[i][j] - a[j][i]));
        if(diff > 1e-1 ) {
          //cout << "WARN: nonsymmetric matrix: " << diff << "\n";
          counter++;
        }
      }
    }
  }
  if( counter> 1)
    return false;
  return true;
}
bool checkskew(complex<double>** a, int n)
{
  //check if the matrix is skew symmetric
  double diff; int counter = 0;
  for(int i=0; i<n; i++) {
    for(int j=0; j<i; j++) {
      if( abs(a[i][j])>1e-6 && abs(a[j][i])>1e-6 ) {
        diff = abs((a[i][j] + a[j][i])/(a[i][j] - a[j][i]));
        if(diff > 2e-2 ) {
          //cout << "WARN: nonsymmetric matrix: " << diff << "\n";
          counter++;
        }
      }
    }
  }
  if( counter> 10)
    return false;
  return true;
}

void copy_m(double** m1, double** m2, int n)
{
  int ds = n*sizeof(double);
  for(int i=0; i<n; i++) memcpy(m1[i], m2[i], ds );
}

void copy_m(complex<double>** m1, complex<double>** m2, int n)
{
  int ds = n*sizeof(complex<double>);
  for(int i=0; i<n; i++) memcpy(m1[i], m2[i], ds );
}

//matrix multiplication(m is multiplied by m1)
void mult(double **m, double **m1, int n)
{
  double **c = createdouble(n);

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      c[i][j] = 0.;

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      for(int k=0; k<n; k++)
        c[i][j] += m1[i][k]*m[k][j];

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      m[i][j] = c[i][j];

  destroy(c,n);
}
void mult(double *m, double **m1, int n)
{
  double *c = new double[n];

  for(int i=0; i<n; i++)
    c[i] = 0.;

  for(int i=0; i<n; i++)
    for(int k=0; k<n; k++)
      c[i] += m1[i][k]*m[k];

  for(int i=0; i<n; i++)
    m[i] = c[i];

  delete[] c;
}

void randomize()
{
  //srand( time(0) );
  long sek;
  time(&sek);
  srand( (unsigned)sek );
  //srand( 0 );
}

//finds the root of a real function in [a,b] using the bisect method
double bisect(double (*f)(double,double,double), double a, double b, double arg1, double arg2)
{
  if( (*f)(a,arg1,arg2)>0. || (*f)(b,arg1,arg2)<0. ) cout << "ERROR: incorrect boundaries\n";

  double am, fam;
  for(int n=0; n<100; n++)
  {
    am = (a+b)/2.;
    fam = (*f)(am, arg1, arg2);
    if( abs(fam)<1e-3 ) break;
    if( fam>0. ) b = am;
    else a = am;
  }
  return am;
}

