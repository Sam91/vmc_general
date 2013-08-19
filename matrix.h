/***************************************************************************
 *   Copyright (C) 2011 by Samuel Bieri   *
 *   samuel.bieri@a3.epfl.ch   *
 *                                                                         *
 ***************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

//compile the Pfaffian functions
#define PFA 0

#include <iostream>
#include <complex>

#define TINY 1e-12

using namespace std;

  //matrix row/col exchange functions
  double update_row      (double*, int, double*, int);
  double update_column   (double*, int, double*, int);
  double det_ratio_row   (double*, int, double*, int);
  double det_ratio_column(double*, int, double*, int);

  complex<double> update_row      (complex<double>*, int, complex<double>*, int);
  complex<double> update_column   (complex<double>*, int, complex<double>*, int);
  complex<double> det_ratio_row   (complex<double>*, int, complex<double>*, int);
  complex<double> det_ratio_column(complex<double>*, int, complex<double>*, int);

  double matrix_element(int, int, int, int, double**);
  complex<double> matrix_element(int, int, int, int, complex<double>**);

  //real matrix functions
  double inverse_m(double**, double**, int);
  int ludcmp(double**, int, int*, double*);
  void lubksb(double**, int, int*, double*);
  double detr_p(double**, int);
  double detr  (double**, int);
  double cofactor(double**, int, int, int);

  //complex matrix functions
  complex<double> inverse_m(complex<double>**, complex<double>**, int);
  int ludcmp(complex<double>**, int, int*, complex<double>*);
  void lubksb(complex<double>**, int, int*, complex<double>*);
  complex<double> detr_p(complex<double>**, int);
  complex<double> detr  (complex<double>**, int);
  complex<double> cofactor(complex<double>**, int, int, int);

#if PFA
  double update_pfa (double*, int, double**, int);
  complex<double> update_pfa (complex<double>*, int, complex<double>**, int);
  double inverse_mpfa(double**, double**, int);
  complex<double> inverse_mpfa(complex<double>**, complex<double>**, int);
#endif

  //matrix inverse functions imported from LAPACK, returning determinant
  double inverse(double*, int);
  complex<double> inverse(complex<double>*, int);
  double inverse(double**, double**, int); //returns det
  complex<double> inverse(complex<double>**, complex<double>**, int); //returns det

  void eigvects(complex<double>**, double*, int);
  void eigvects(double**, double*, int);
  void eigvals(complex<double>**, double*, int);
  void eigvals(double**, double*, int);

  double inverse_mkeep(double**, double**, int);
  complex<double> inverse_mkeep(complex<double>**, complex<double>**, int);
#endif

