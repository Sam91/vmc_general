/***************************************************************************
 *   Copyright (C) 2011 by Samuel Bieri                                    *
 *   samuel.bieri@a3.epfl.ch                                               *
 *                                                                         *
 ***************************************************************************/

/* This source implements all matrix related methods of matrix.h */

#include "matrix.h"
#include "helperfunctions.h"

#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>

//Wimmer's Pfaffian package
#if PFA
#include "pfapack.h"
#endif

/* Calculates the inverse of a and returns the det. Destroys the matrix a. The inv is written to y */
double inverse_m(double **a, double **y, int n)
{
  int i,j, *indx;
  double d;
  indx = new int[n];
  double *col = new double[n];

  //cout << "inverse_m (real)\n";

  if( ludcmp(a,n,indx,&d)!=0 )
    d = 0.;
  else {
    //write_m(a,n);
    for (j=0; j<n; j++) {
      d *= a[j][j];
      for (i=0; i<n; i++) col[i]=0.;
      col[j]=1.;
      lubksb(a,n,indx,col);
      for (i=0; i<n; i++) y[i][j]=col[i];
    }
  }
  delete[] indx;
  delete[] col;
  return d;
}

complex<double> inverse_m(complex<double> **a, complex<double> **y, int n)
{
  int i,j, *indx;
  complex<double> d;
  indx = new int[n];
  complex<double> *col = new complex<double>[n];

  //cout << "inverse_m (complex)\n";

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      y[i][j] = 0.;

  if( ludcmp(a,n,indx,&d)!=0 )
    d = 0.;
  else {
    for(j=0;j<n;j++) {
      d *= a[j][j];
      for(i=0;i<n;i++) col[i]=complex<double>(0.,0.);
      col[j]=complex<double>(1.,0.);
      lubksb(a,n,indx,col);
      for(i=0;i<n;i++) y[i][j]=col[i];
    }
  }
  delete[] indx;
  delete[] col;
  return d;
}

/* Calculates the inverse of a skew symmetric matrix a and returns the Pfaffian. a is destroyed and y returns the inv*/
#if PFA
double inverse_mpfa(double **a, double **y, int n)
{
  int i,j,info=0, LWORK;

  if( n<2 || n%2==1 ) {
    cout << "WARN: inapropriate matrix size for Pfaffian calculation\n";
    return 0.;
  }

  vector<int> IWORK(n+1);
  vector<double> WORK(n+1);
  vector<double> A(n*n);

  //fill the L part of the A with values in fortran convention
  fill(A.begin(), A.end(), 0.);
  for(i=0;i<n;i++) {
    for(j=0;j<i;j++) {
      A[i+j*n] = a[i][j];
    }
  }

  //Do the workspace quer first
  LWORK = -1;
  dsktrf_("L", "P", &n, &A[0], &n, &IWORK[0], &WORK[0], &LWORK, &info);

  if(info) {
    cout<< "WARN: There is a problem in the workspace call\n";
    return 0.;
  }

  LWORK=static_cast<int>(WORK[0]);
  if(LWORK>n+1) WORK.resize(LWORK); 

  //compute the partial tri-diagonal form of the matrix
  dsktrf_("L", "P", &n, &A[0], &n, &IWORK[0], &WORK[0], &LWORK, &info);

  //accumulate the Pfaffian
  double pfa = 1.;
  for(i=1; i<n; i+=2)
  {
    pfa *= -A[i+(i-1)*n];
    if( IWORK[i] != i+1 ) pfa *= -1.;
  }

  double apfa = abs(pfa);
  if( apfa < 1e-4 ) return 0.;
  //if( apfa==std::numeric_limits<double>::infinity() ) {cout << "Warn: determinant overflow (" << pfa << "); returning zero\n"; return 0.;}
  if( apfa > 1e100 ) {cout << "Warn: determinant overflow (" << std::scientific << pfa << "); returning zero\n"; return 0.;}

//  double *A1 = new double[N*N];
//  double **bkpa = createdouble(n);
//  copy_m(bkpa, a, n);

  //Next, we compute the inverse matrix from scratch (it would be slightly more efficient to use the LTL decomposition, but I don't do this for the moment)
  //double det = inverse_m(a, y, n);
  inverse_m(a, y, n);


//  if( !checkskew(y, n) )
//  {
//    cout << "Non-skew! pfa= " << pfa << "; det=" << det << "\n";
//    cout << "matrix:\n";
//    write_m(bkpa, n);
//    save_m(bkpa,n,"d3");
//    cout << "inv:\n";
//    write_m(y, n);
//  } else
//    cout << "Skew ok after inverse_mpfa\n";
//  destroy(bkpa,n);

  //anti-symmetrize
//  for(int i=0; i<n; i++) {
//    y[i][i] = 0.;
//    for(int j=0; j<i; j++) {
//      y[i][j] = (y[i][j]-y[j][i])/2.;
//      y[j][i] = -y[i][j];
//    }
//  }

  return pfa;
}

complex<double> inverse_mpfa(complex<double> **a, complex<double> **y, int n)
{
  complex<double> pfa;
  int i,j,info=0,LWORK;

  if( n<2 || n%2==1 ){
    cout << "WARN: inapropriate matrix size for Pfaffian calculation: "<< n << "\n";
    return 0.;
  }

  vector<int> IWORK(n+1);
  vector<complex<double> > WORK(n+1);
  vector<complex<double> > A(n*n);

  //fill the L part of the A with values in fortran convention
  fill(A.begin(), A.end(), 0.0);
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      A[i+j*n] = a[i][j];

  //Do the workspace quer first
  LWORK = -1;
  zsktrf_("L", "P", &n, &A[0], &n, &IWORK[0], &WORK[0], &LWORK, &info);

  if(info) {
    cout<< "WARN: There is a problem in the workspace call\n";
    return 0.;
  }

  LWORK=static_cast<int>(real(WORK[0]));
  if(LWORK>n+1) WORK.resize(LWORK);

  //compute the partial tri-diagonal form of the matrix
  zsktrf_("L", "P", &n, &A[0], &n, &IWORK[0], &WORK[0], &LWORK, &info);

  //accumulate the Pfaffian
  pfa = 1.;
  for(i=1; i<n; i+=2)
  {
    pfa *= -A[i+(i-1)*n];
    if( IWORK[i] != i+1 ) pfa *= -1.;
  }

  double apfa = abs(pfa);
  if( apfa < 1e-4 ) return 0.;
  if( apfa > 1e100 ) {cout << "Warn: determinant overflow (" << std::scientific << pfa << "); returning zero\n"; return 0.;}

  //Next, we compute the inverse matrix from scratch (it would be slightly more efficient to use the LTL decomposition, but I don't do this for the moment)
  inverse_m(a, y, n);

  return pfa;
}
#endif

int ludcmp(double **a, int n, int *indx, double *d)
{
  int i,j,k;
  int imax = -1;
  double big,dum,sum,temp;
  double *vv;

  //cout << "ludcmp (real)\n";

  vv = new double[n];
  (*d)=1.;
  for(i=0; i<n; i++) {
    big=0.;
    for(j=0; j<n; j++)
      if( (temp=fabs(a[i][j]))>big ) big = temp;
    if( big == 0. ) {
      cout << "Singular matrix in routine ludcmp (real)\n"; //exit(-1);
      (*d)=0.; return -1;
    }
    vv[i]=1./big;
  }

  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      (*d) = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if( a[j][j] == 0.) a[j][j] = TINY;
    if (j != n-1) {
      dum=1./(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  delete[] vv;
  return 0;
}

int ludcmp(complex<double> **a, int n, int *indx, complex<double> *d)
{
  int i,j,k;
  int imax = -1;
  double dd, temp, big, dum;
  complex<double> dum2, sum;
  double *vv;

  //cout << "ludcmp (complex)\n";

  vv = new double[n];
  dd = 1.;
  for (i=0;i<n;i++) {
    big = 0.;
    for(j=0; j<n; j++)
      if ((temp=abs(a[i][j])) > big) big=temp;
    if( big == 0. ) {
      cout << "Singular matrix in routine ludcmp (complex)\n"; //exit(-1);
      (*d)=0.; delete[] vv; return -1;
    }
    vv[i]=1./big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for( k=0;k<i;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if( (dum=vv[i]*abs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum2=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum2;
      }
      dd = -dd;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    //Not sure why we need the following code
    if (a[j][j] == 0.) a[j][j] = TINY;
    if (j != n-1) {
      dum2 = pow(a[j][j],-1);
      for (i=j+1;i<n;i++) a[i][j] *= dum2;
    }
  }
  delete[] vv;
  (*d) = complex<double>(dd,0.);
  return 0;
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=-1,ip,j;
  double sum;

  //cout << "lubksb (real)\n";

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if( ii>=0 )
      for (j=ii; j<i; j++)
        sum -= a[i][j]*b[j];
    else if( sum ) ii=i;
    b[i]=sum;
  }
  for(i=n-1; i>=0; i--) {
    sum=b[i];
    for(j=i+1; j<n; j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void lubksb(complex<double> **a, int n, int *indx, complex<double> b[])
{
  int i,ii=-1,ip,j;
  complex<double> sum;

  //cout << "lubksb (complex)\n";

  for( i=0;i<n;i++ ) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if( ii>=0 )
      for (j=ii; j<i; j++)
        sum -= a[i][j]*b[j];
    else if( sum!= 0. ) ii=i;
    b[i] = sum;
  }
  for( i=n-1; i>=0; i-- ) {
    sum=b[i];
    for( j=i+1; j<n; j++ ) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

/* Calculates the det. Keeps the original matrix a. */
double detr_p(double **a, int n)
{
  int i;
  double  ret;
  double **a2;

  a2 = new ( double(*[n]) ) ;
  for(i=0; i<n; i++) {
    a2[i] = new double[n];
    for(int j=0; j<n; j++) a2[i][j] = a[i][j];
  }
  ret = detr(a2,n);

  for(i=0;i<n;i++) delete[] a2[i]; delete[] a2;

  return ret;
}
complex<double> detr_p(complex<double> **a, int n)
{
  int i;
  complex<double>  ret;
  complex<double> **a2;

  a2 = new ( complex<double>(*[n]) ) ;
  for(i=0; i<n; i++) {
    a2[i] = new complex<double>[n];
    for(int j=0; j<n; j++) a2[i][j] = a[i][j];
  }
  ret = detr(a2,n);

  for(i=0;i<n;i++) delete[] a2[i]; delete[] a2;

  return ret;
}

double detr(double** a, int n)
{
  int j,*indx;
  double d;
  indx = new int[n];

  ludcmp(a,n,indx,&d);

  for (j=0;j<n;j++) d *= a[j][j];

  delete[] indx;
  return d;
}
complex<double> detr(complex<double>** a, int n)
{
  int j,*indx;
  complex<double> d;
  indx = new int[n];

  ludcmp(a,n,indx,&d);

  for( j=0;j<n;j++ ) d *= a[j][j];

  delete[] indx;
  return d;
}

double cofactor(double** a, int n, int i, int j)
{
//  cout << "Entering cofactor matrix with n=" << n << ", i="<< i << ", j="<< j << "\n";

  double **subm;
  subm = new(double(*[n]) );
  int ii,jj;

  ii=0;
  for(int k=0; k<n; k++) {
    if(k==i) continue;
    subm[ii] = new double[n];
    jj=1;
    for(int l=0;l<n;l++) if(l!=j) subm[ii][jj++] = a[k][l];
    ii++;
  }
  //write_m(subm,n);

  double det = detr(subm,n);

  for(int k=0; k<n; k++) delete[] subm[k]; delete[] subm;

  return (((i+j)%2)==0)?det:-det;
}

complex<double> cofactor(complex<double>** a, int n, int i, int j)
{
  complex<double> **subm;
  subm = new(complex<double>(*[n]) );
  int ii,jj;

  ii=0;
  for(int k=0; k<n; k++) {
    if(k==i) continue;
    subm[ii] = new complex<double>[n];
    jj=1;
    for(int l=0;l<n;l++) if(l!=j) subm[ii][jj++] = a[k][l];
    ii++;
  }

  complex<double> det = detr(subm,n);

  for(int k=0; k<n; k++) delete[] subm[k]; delete[] subm;

  return (((i+j)%2)==0)?det:-det;
}

/* Sherman-Morrison algorithm to compute the inverse of matrix with one row/col replaced */
/* See "Numerical Recipes" for more details                                                */
double update_row( double *row, int row_number, double *inv, int n)
{
  /*Calculate the ratio of the determinants */
  double r = det_ratio_row(row, row_number, inv, n);

  /* checking for degeneracy */
  if ( fabs(r)<TINY ) {
//    ZERO_STATUS = 1; cout << "ZS by row\n";
    return( 0. ); //We return 0 here and the move will be rejected
  }

  double x;

  for(int j=0; j<n; j++)
  {
    if( j!=row_number ) {
      x = 0.;
      for(int i=0; i<n; i++) x += row[i]*inv[i*n+j];
      x/=r;
      for(int i=0; i<n; i++) inv[i*n+j] -= inv[i*n+row_number]*x;
    }
  }
  for(int i=0; i<n; i++) inv[i*n+row_number]/=r;
  return r;
}

complex<double> update_row( complex<double> *row, int row_number, complex<double> *inv, int n)
{
  complex<double> r = det_ratio_row(row, row_number, inv, n);

  /* checking for degeneracy */
  if ( abs(r)<TINY ) {
//    ZERO_STATUS = 1; //cout << "ZS by row\n";
    return( 0. ); //We return 0 here and the move will be rejected
  }

  complex<double> x;

  for(int j=0; j<n; j++)
  {
    if (j!=row_number) {
      x = 0.;
      for(int i=0; i<n; i++) x += row[i]*inv[i*n+j];
      x/=r;
      for(int i=0; i<n; i++) inv[i*n+j] -= inv[i*n+row_number]*x;
    }
  }
  for(int i=0; i<n; i++) inv[i*n+row_number]/=r;
  return r;
}

double update_column( double *column, int column_number, double *inv, int n)
{
  double r = det_ratio_column(column, column_number, inv, n);

  /* checking for degeneracy */
  if( fabs(r)<TINY ) {
//    ZERO_STATUS = 1; //cout << "ZS by col\n";
    return( 0. ); //We return 0 here and the det will be calculated from scratch.
  }

  double x;

  for(int j=0; j<n; j++)
  {
    if( j!=column_number ) {
      x = 0.;
      for(int i=0; i<n; i++) x += column[i]*inv[j*n+i];
      x/=r;
      for(int i=0; i<n; i++) inv[j*n+i] -= inv[column_number*n+i]*x;
    }
  }
  for(int i=0; i<n; i++) inv[column_number*n+i]/=r;
  return r;
}

complex<double> update_column( complex<double> *column, int column_number, complex<double> *inv, int n)
{
//  cout << "Update column\n";
//check this code
/*  complex<double> **dm = createcomplex(n);
  inverse_mkeep(inv, dm, n);
  for(int i=0; i<n; i++) dm[i][column_number] = column[i];
  complex<double> **dm2 = createcomplex(n);
*/
  complex<double> r = det_ratio_column(column, column_number, inv, n);

  /* checking for degeneracy */
  if( abs(r)<TINY ) {
    return( 0. ); //We return 0 here and the det will be calculated from scratch.
  }

  complex<double> x;

  for(int j=0; j<n; j++)
  {
    if (j!=column_number) {
      x = 0.;
      //for(i=0; i<n; i++) x += column[i]*inv[j][i];
      for(int i=0; i<n; i++) x += column[i]*inv[j*n+i];
      x/=r;
      //for(i=0; i<n; i++) inv[j][i] -= inv[column_number][i]*x;
      for(int i=0; i<n; i++) inv[j*n+i] -= inv[column_number*n+i]*x;
    }
  }
  //for(i=0; i<n; i++) inv[column_number][i]/=r;
  for(int i=0; i<n; i++) inv[column_number*n+i]/=r;

/*  inverse_mkeep(inv, dm2, n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if( abs(dm2[i][j]-dm[i][j])>TINY ) cout << "ERROR in update_column\n";
  destroy(dm,n); destroy(dm2,n);
*/
  return r;
}

/* Updates one row and column of the inverse of a skew symmetric matrix */
#if PFA
double update_pfa( double *column, int col_number, double **inv, int n)
{
  int i,j;
  double x;

  //Pf(B)/Pf(A) = det(\tilde A)/det(A), where \tilde A has only one column replaced

/*  double **a1 = createdouble(n);
  double **a2 = createdouble(n);
  copy_m(a1, inv, n);
*/
/*  if( !checkskew(inv, n) ) {
    cout << "Skew matrix update problem before start...\n";
    write_m(inv, n);
  }
*/

  double r = 0.;
  for (i=0; i<n; i++) r += column[i]*inv[col_number][i];

  double fr = abs(r);
  if ( fr<1e-3 )
    return 0.;

//  if( fr==std::numeric_limits<double>::infinity() || fr > 1e100) {
  if( fr > 1e100) {
    cout << "update_pfa: determinant overflow (real)\n";
    return 0.; //We return 0 here and the move will be rejected
  }

  //do the column replacement
  for (j=0; j<n; j++) {
    if (j!=col_number) {
      x = 0.;
      for (i=0;i<n;i++) x += column[i]*inv[j][i];
      x /= r;
      for (i=0;i<n;i++) inv[j][i] -= inv[col_number][i]*x;
    }
  }
  for (i=0; i<n; i++) inv[col_number][i]/=r;

  //copy_m(a2, inv, n);

  //double r2 = r;
  //double r2 = 0.;
  //for (i=0; i<n; i++) r2 -= column[i]*inv[i][col_number];

  //do the row replacement
  for (j=0; j<n; j++) {
    if (j!=col_number) {
      x = 0.;
      for (i=0; i<n; i++) x -= column[i]*inv[i][j];
      x /= r;
      for (i=0; i<n; i++) inv[i][j] -= inv[i][col_number]*x;
    }
  }
  for (i=0; i<n; i++) inv[i][col_number]/=r;


/*  if( !checkskew(inv, n) ) {
    cout << "Skew matrix update problem. r1=" << r << "; r2=" << r2 << "\n";
    cout << "original:\n";
    write_m(a1, n);
    save_m(a1, n, "inv1");

    cout << "intermediate inv:\n";
    write_m(a2,n);

    cout << "new:\n";
    write_m(inv, n);
    save_m(inv, n, "inv2");

    cout << "replacement vector (" << col_number << ": {";
    for(i=0; i<n; i++) cout << column[i] << ", ";
    cout << "}\n";
    save_v(column, n, "col");
  }
*/

//  if( abs(r/r2-1.)>1e-3 ) 
//    cout << "r1, r2: " << r << ", " << r2 << "\n";

  //anti-symmetrize
  for(i=0; i<n; i++) {
    inv[i][i] = 0.;
    for(j=0; j<i; j++) {
      inv[i][j] = (inv[i][j]-inv[j][i])/2.;
      inv[j][i] = -inv[i][j];
    }
  }

  //destroy(a1, n);
  //destroy(a2, n);
  return r;
}

complex<double> update_pfa( complex<double> *column, int col_number, complex<double> **inv, int n )
{
  int i,j;

//  complex<double> **a = createcomplex(n);
//  copy_m(a, inv, n);

  //Pf(B)/Pf(A) = det(\tilde A)/det(A), where \tilde A has only one column replaced
  complex<double> r = 0.;
  for (i=0; i<n; i++) r += column[i]*inv[col_number][i];

  double fr = abs(r);
  if ( fr<1e-3 )
    return 0.;

//  if( fr==std::numeric_limits<double>::infinity() || fr > 1e100) {
  if( fr > 1e100 ) {
    cout << "update_pfa: determinant overflow (complex)\n";
    return 0.; //We return 0 here and the move will be rejected
  }

  //do the column replacement
  complex<double> x;
  for( j=0;j<n;j++ ) {
    if( j!=col_number ) {
      x = 0.;
      for( i=0;i<n;i++ ) x += column[i]*inv[j][i];
      x/=r;
      for( i=0;i<n;i++ ) inv[j][i] -= inv[col_number][i]*x;
    }
  }
  for (i=0;i<n;i++) inv[col_number][i]/=r;

  //do the row replacement
  //complex<double> r2 = r;
  //for (i=0; i<n; i++) r2 -= column[i]*inv[i][col_number];

  for (j=0; j<n; j++) {
    if (j!=col_number) {
      x = 0.;
      for (i=0; i<n; i++) x -= column[i]*inv[i][j];
      x/=r;
      for (i=0; i<n; i++) inv[i][j] -= inv[i][col_number]*x;
    }
  }
  for (i=0; i<n; i++) inv[i][col_number]/=r;

/*  if( !checkskew(inv, n) ) {
    cout << "Skew matrix update problem\n";
    cout << "original:\n";
    write_m(a, n);
    //save_m(a, n, "inv1");

    cout << "new:\n";
    write_m(inv, n);
    //save_m(inv, n, "inv2");

    cout << "replacement vector (" << col_number << ": {";
    for(i=0; i<n; i++) cout << column[i] << ", ";
    cout << "}\n";
    //save_v(column, n, "col");
  }
  destroy(a, n);
*/

  //anti-symmetrize
  for(i=0; i<n; i++) {
    inv[i][i] = 0.;
    for(j=0; j<i; j++) {
      inv[i][j] = (inv[i][j]-inv[j][i])/2.;
      inv[j][i] = -inv[i][j];
    }
  }

  //return r here; it is the ratio of the pfaffians
  return r;
}
#endif

double det_ratio_row( double *row, int row_number, double *inv, int n)
{
  double r = 0.;
  for(int i=0; i<n; i++) r += row[i]*inv[i*n+row_number];
  return r;
}

complex<double> det_ratio_row( complex<double> *row, int row_number, complex<double> *inv, int n)
{
  complex<double> r = 0.;
  for(int i=0; i<n; i++) r += row[i]*inv[i*n+row_number];
  return r;
}

double det_ratio_column( double *column, int column_number, double *inv, int n)
{
  double r = 0.;
  for(int i=0; i<n; i++) r += column[i]*inv[column_number*n+i];
  return r;
}

complex<double> det_ratio_column( complex<double> *column, int column_number, complex<double> *inv, int n)
{
  complex<double> r = 0.; //int cn = column_number*n;
  for(int i=0; i<n; i++) r += column[i]*inv[column_number*n+i];
  return r;
}


/* LAPACK inverse function */

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    void zgetrf_(int* M, int *N, complex<double>* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
    void zgetri_(int* N, complex<double>* A, int* lda, int* IPIV, complex<double>* WORK, int* lwork, int* INFO);

/* ZHEEV prototype */
    void zheev_( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info );
    void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
}

double inverse(double* A, int n)
{
  int *IPIV = new int[n+1];
  int LWORK = n*n;
  double *WORK = new double[LWORK];
  int INFO;

  dgetrf_(&n,&n,A,&n,IPIV,&INFO);

  if(INFO > 0) { //singular matrix
    delete[] IPIV; delete[] WORK;
    return 0.;
  }

  double det = 1.; bool sgn = true;
  for(int i=0; i<n; i++) {
    det *= A[(1+n)*i];
    if( IPIV[i] != (i+1) ) sgn = !sgn;
  }

  dgetri_(&n,A,&n,IPIV,WORK,&LWORK,&INFO);

  delete[] IPIV; delete[] WORK;

  return sgn?det:-det;
}

double inverse(double** mat, double** inv, int n)
{
  int sd = n*sizeof(double);
  double *A = new double[n*n];
  for(int i=0; i<n; i++) memcpy(A + n*i, mat[i], sd );

  int *IPIV = new int[n+1];
  int LWORK = n*n;
  double *WORK = new double[LWORK];
  int INFO;

  dgetrf_(&n, &n, A, &n, IPIV, &INFO);

  if(INFO > 0) { //singular matrix
    delete[] A; delete[] IPIV; delete[] WORK;
    return 0.;
  }

  double det = 1.; bool sgn = true;
  for(int i=0; i<n; i++) {
    det *= A[(1+n)*i];
    if( IPIV[i] != (i+1) ) sgn = !sgn;
  }
  
  dgetri_(&n, A, &n, IPIV, WORK, &LWORK, &INFO);

  for(int i=0; i<n; i++) memcpy(inv[i], A + n*i, sd );

  delete[] A; delete[] IPIV; delete[] WORK;

  return sgn?det:-det;
}

complex<double> inverse(complex<double>* A, int N)
{
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  complex<double> *WORK = new complex<double>[LWORK];
  int INFO;

  zgetrf_(&N,&N,A,&N,IPIV,&INFO);

  if(INFO > 0) { //singular matrix
    delete[] IPIV; delete[] WORK;
    return 0.;
  }

  complex<double> det = 1.; bool sgn = true;
  for(int i=0; i<N; i++) {
    det *= A[(1+N)*i];
    if( IPIV[i] != (i+1) ) sgn = !sgn;
  }

  zgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

  delete[] IPIV; delete[] WORK;

  return sgn?det:-det;
}

complex<double> inverse(complex<double>** mat, complex<double>** inv, int N)
{
  int sd = N*sizeof(complex<double>);
  complex<double> *A = new complex<double>[N*N];
  for(int i=0; i<N; i++) memcpy(A + N*i, mat[i], sd );

  int *IPIV = new int[N+1];
  int LWORK = N*N;
  complex<double> *WORK = new complex<double>[LWORK];
  int INFO;

  zgetrf_(&N, &N, A, &N, IPIV, &INFO);

  if(INFO > 0) { //singular matrix
    delete[] A; delete[] IPIV; delete[] WORK;
    return 0.;
  }

  complex<double> det = 1.; bool sgn = true;
  for(int i=0; i<N; i++) {
    det *= A[(1+N)*i];
    if( IPIV[i] != (i+1) ) sgn = !sgn;
  }

  zgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

  for(int i=0; i<N; i++) memcpy(inv[i], A + N*i, sd );

  delete[] A; delete[] IPIV; delete[] WORK;

  return sgn?det:-det;
}

//calculate the eigenvalues and eigenvectors of a complex hermitian matrix
//eigenvectors are returned as rows in m, i.e., m[j][n] where n is the eigenvalue index; they are orthogonal
void eigvects(complex<double> **m, double* w, int N)
{
  int lwork, info, n, lda;

  //cout << "eigenvects with:\n";
  //write_m(m, N);

  dcomplex wkopt;
  dcomplex *work;
  double *e = new double[N];
  double *rwork = new double[3*N-2];

  dcomplex *AT = new dcomplex[N*N];

  for(int i=0; i<N; i++) {
    for(int j=0; j<=i; j++) {
      AT[i+N*j].re = m[i][j].real();
      AT[i+N*j].im = m[i][j].imag();
    }
  }

  n = N; lda = N;
  lwork = -1;

  zheev_(const_cast<char *>("V"), const_cast<char *>("L"), &n, AT, &lda, w, &wkopt, &lwork, rwork, &info);

  lwork = (int)wkopt.re;
  work = new dcomplex[lwork]; 

  zheev_(const_cast<char *>("V"), const_cast<char *>("L"), &n, AT, &lda, e, work, &lwork, rwork, &info);

  if( info>0 ) cout << "ERROR in zheev routine\n";

  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      m[i][j] = complex<double>(AT[j+N*i].re, AT[j+N*i].im);
    }
//    cout << "\n";
  }

  for(int i=0; i<N; i++) w[i]=e[i];

  delete[] work; delete[] AT; delete[] e; delete[] rwork;
}

void eigvects(double **m, double* w, int N)
{
  int lwork, info, n, lda;

  //cout << "eigenvects with:\n";
  //write_m(m, N);

  double wkopt;
  double *work;

  double *AT = new double[N*N];

  for(int i=0; i<N; i++) {
    for(int j=0; j<=i; j++) {
      AT[i+N*j] = m[i][j];
    }
  }

  n = N; lda = N;
  lwork = -1;

  dsyev_(const_cast<char *>("V"), const_cast<char *>("L"), &n, AT, &lda, w, &wkopt, &lwork, &info);

  lwork = (int)wkopt;
  work = new double[lwork];

  dsyev_(const_cast<char *>("V"), const_cast<char *>("L"), &n, AT, &lda, w, work, &lwork, &info);

  if( info>0 ) cout << "ERROR in dsyev routine\n";

  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      m[i][j] = AT[j+N*i];
    }
  }

  delete[] work; delete[] AT;
}

void eigvals(complex<double> **m, double* w, int N)
{
  int lwork, info, n, lda;

  //cout << "eigenvals(c) with:\n";
  //write_m(m, N);

  dcomplex wkopt;
  dcomplex *work;
  double *rwork = new double[3*N-2];

  dcomplex *AT = new dcomplex[N*N];
  //cout << "starting copy loop\n";

  for(int i=0; i<N; i++) {
    //cout << i << "\n"; 
    for(int j=0; j<=i; j++) {
      //cout << j << "; "; 
      AT[i+N*j].re = m[i][j].real();
      AT[i+N*j].im = m[i][j].imag();
    }
  }

  n = N; lda = N;
  lwork = -1;

  zheev_(const_cast<char *>("N"), const_cast<char *>("L"), &n, AT, &lda, w, &wkopt, &lwork, rwork, &info);

  lwork = (int)wkopt.re;
  work = new dcomplex[lwork];

  zheev_(const_cast<char *>("N"), const_cast<char *>("L"), &n, AT, &lda, w, work, &lwork, rwork, &info);

  if( info>0 ) cout << "ERROR in zheev routine\n";

  delete[] work; delete[] AT; delete[] rwork;
}

void eigvals(double **m, double* w, int N)
{
  int lwork, info, n, lda;

  //cout << "eigenvals(d) with:\n";
  //write_m(m, N);

  double wkopt;
  double *work;

  double *AT = new double[N*N];

  for(int i=0; i<N; i++) {
    for(int j=0; j<=i; j++) {
      AT[i+N*j] = m[i][j];
    }
  }

  n = N; lda = N;
  lwork = -1;

  dsyev_(const_cast<char *>("N"), const_cast<char *>("L"), &n, AT, &lda, w, &wkopt, &lwork, &info);

  lwork = (int)wkopt;
  work = new double[lwork];

  dsyev_(const_cast<char *>("N"), const_cast<char *>("L"), &n, AT, &lda, w, work, &lwork, &info);

  if( info>0 ) cout << "ERROR in dsyev_ routine\n";

  delete[] work; delete[] AT;
}

//inverse calculation that keeps the original matrix
double inverse_mkeep(double **a, double **y, int n)
{
  double **m = createdouble(n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      m[i][j] = a[i][j];
  double r = inverse_m(m,y,n);
  destroy(m,n);
  return r;
}
complex<double> inverse_mkeep(complex<double> **a, complex<double> **y, int n)
{
  complex<double> **m = createcomplex(n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      m[i][j] = a[i][j];
  complex<double> r = inverse_m(m,y,n);
  destroy(m,n);
  return r;
}


//some test function for inversion and diagonalization speed
void test_diag( int n, int N )
{
  cout << "Starting diagonalization test (n=" << n << ")...\n";

  double t0, dt1, dt2, dt3, dt4;
  double* v = new double[n];

  double** mr = createdouble( n );
  double** mr1 = createdouble( n );
  double** mr2 = createdouble( n );
  double** mr3 = createdouble( n );
  double* r = new double[n*n];
  double detr1, detr2;

  dt1=0.; dt2=0.; dt3=0.; dt4=0.;
  for(int i=0; i<N; i++)
  {
    for(int k1=0; k1<n; k1++)
      for(int k2=0; k2<n; k2++) {
        mr[k1][k2] = 2.*GET_RAND;
        //mr[k2][k1] = mr[k1][k2];
        r[k1+n*k2] = mr[k1][k2];
        //r[k2+n*k1] = mr[k2][k1];
      }

    t0 = time(0);
    detr2 = inverse(mr, mr2, n);
    //inverse(r, n);
    dt1 += time(0) - t0;

    copy_m(mr1, mr, n);
    t0 = time(0);
    detr1 = inverse_m(mr1, mr3, n);
    dt2 += time(0) - t0;

    if( abs(detr1/detr2-1.)>1e-10 ) cout << "WARN: det mismatch: " << detr1 << ", " << detr2 << "\n";

    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        if( abs(mr2[i][j]-mr3[i][j])>1e-3 ) cout << "Warn: r error\n";

    copy_m(mr1, mr, n);
    t0 = time(0);
    eigvals(mr1, v, n);
    dt3 += time(0) - t0;

    copy_m(mr1, mr, n);
    t0 = time(0);
    eigvects(mr1, v, n);
    dt4 += time(0) - t0;
  }

  cout << "Lapack inversion of real matrix: " << std::fixed << setprecision(3) << dt1/(double)N << " sec.\n";
  cout << "In-house inversion of real matrix: " << std::fixed << setprecision(3) << dt2/(double)N << " sec.\n";
  cout << "Lapack full diag of real matrix: " << std::fixed << setprecision(3) << dt3/(double)N << " sec.\n";
  cout << "Lapack eigenvects of real matrix: " << std::fixed << setprecision(3) << dt4/(double)N << " sec.\n";

  destroy( mr, n);
  destroy( mr1, n);
  destroy( mr2, n);
  destroy( mr3, n);
  delete[] r;

  cout << "\n";

  complex<double>** mc = createcomplex( n );
  complex<double>** mc1 = createcomplex( n );
  complex<double>** mc2 = createcomplex( n );
  complex<double>** mc3 = createcomplex( n );
  complex<double>* c = new complex<double>[n*n];
  complex<double> detc1, detc2;

  dt1=0.; dt2=0.; dt3=0.; dt4=0.;
  for(int i=0; i<N; i++)
  {
    //cout << "Creating matrix...";
    for(int k1=0; k1<n; k1++) {
      for(int k2=0; k2<n; k2++) {
        mc[k1][k2] = 2.*complex<double>(GET_RAND, GET_RAND);
        //mc[k2][k1] = conj(mc[k1][k2]);
        c[k1+n*k2] = mc[k1][k2];
        //c[k2+n*k1] = mc[k2][k1];
      }
    }
    //cout << "done.\n";
    
    t0 = time(0);
    detc1 = inverse(mc, mc2, n);
    //inverse(c, n);
    dt1 += time(0) - t0;

    copy_m(mc1, mc, n);
    t0 = time(0);
    detc2 = inverse_m(mc1, mc3, n);
    dt2 += time(0) - t0;

    if( abs(detc1/detc2-1.)>1e-10 ) cout << "WARN: det mismatch: " << detc1 << ", " << detc2 << "\n";

    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        if( abs(mc2[i][j]-mc3[i][j])>1e-3 ) cout << "Warn: c error\n";

    copy_m(mc1, mc, n);
    t0 = time(0);
    eigvals(mc1, v, n);
    dt3 += time(0) - t0;

    copy_m(mc1, mc, n);
    t0 = time(0);
    eigvects(mc1, v, n);
    dt4 += time(0) - t0;
  }

  cout << "Lapack inversion of complex matrix: " << std::fixed << setprecision(3) << dt1/(double)N << " sec.\n";
  cout << "In-house inversion of complex matrix: " << std::fixed << setprecision(3) << dt2/(double)N << " sec.\n";
  cout << "Lapack full diag of complex matrix: " << std::fixed << setprecision(3) << dt3/(double)N << " sec.\n";
  cout << "Lapack eigenvects of complex matrix: " << std::fixed << setprecision(3) << dt4/(double)N << " sec.\n";

  destroy( mc, n);
  destroy( mc1, n);
  destroy( mc2, n);
  destroy( mc3, n);
  delete[] c;

  delete[] v;
}

