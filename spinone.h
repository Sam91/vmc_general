#ifndef SPINONE_H
#define SPINONE_H

#include "helperfunctions.h"

//Helper class containing data which specify a spin-one state (single site) as well as some usefull methods
class spinone
{

public:
  spinone(); //creates a spin one state (0,0,1),(0,1,0)
  spinone(double* u, double* v); //generic u and v vectors
  spinone(double eta, double phi, double theta, double psi); //spin one state corresponding to rotations by Euler angles
  spinone(const spinone&);
  void operator = (const spinone&);

  ~spinone();

  void init();

  //calculate the (3x3) D matrix corresponding to a given spin one state
  void createD();

  //rotate the current u and v vectors by Euler angles
  void rotate(double phi, double theta, double psi);

  //calculate the ground and highest eigenstates of the D operator in Fermionic space
  void projectD();

  void print(); //print our the data
  void normalize(complex<double> *vect, int n); //ortho-normalize the u- and v-vectors

  complex<double> **D; //D-matrix resulting from the u- and v-vectors

  double getU(int);
  double getV(int);

  //does the D matrix corresponding to this spin state couple the z-fermion or not
  bool uncoupled_z();

private:

  //real and imaginary part of the director
  double *u;
  double *v;

};

#endif
