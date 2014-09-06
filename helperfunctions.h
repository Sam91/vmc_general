/***************************************************************************
 *   Copyright (C) 2011 by Samuel Bieri   *
 *   samuel.bieri@a3.epfl.ch   *
 *                                                                         *
 ***************************************************************************/

#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <complex>
#include "randoma.h"

using namespace std;

#include <stdlib.h>
#ifndef GET_RAND
#define GET_RAND ( (double)rand()/((double)RAND_MAX + 1.0) ) // the standard c library in stdlib.h
#endif


//2D array creation functions
int **createint(int);
int **createint(int,int);
double **createdouble(int);
double **createdouble(int,int);
complex<double> **createcomplex(int);
complex<double> **createcomplex(int,int);
complex<double> ***createcomplex(int,int,int);
void destroy(int**, int);
void destroy(double**, int);
void destroy(double***, int, int);
void destroy(complex<double>**, int);

//matrix print functions
void write_m(double**, int);
void write_m(complex<double>**, int);
void write_m(int**, int);
void write_m(double*, int);
void write_m(complex<double>*, int);

void save_m(double**, int, const char*);
void load_m(double**, int, const char*);
int save_v(double*, int, const char*);
int load_v(double*, int, const char*);
int save_v(int*, int, const char*);
int load_v(int*, int, const char*);

//int fappend(double*, int, const char*);
int fappend(const double*, int, string);

void copy_m(double**, double**, int );
void copy_m(complex<double>**, complex<double>**, int );

void estimate_time( double time_start, int to_do, int done );

bool checkskew(double**, int);
bool checkskew(complex<double>**, int);

void mult(double**, double**, int n);
void mult(double*, double**, int n);

void randomize();

void test_diag(int,int);
double bisect(double (*f)(double,double,double), double a, double b, double arg1, double arg2);

int factorial(int);

#endif

