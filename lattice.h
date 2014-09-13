//define a general D-dimensional bravais lattice

#ifndef LATTICE_H
#define LATTICE_H

//dimension of the lattice
#define DIM 2

//number of sublattice sites
#define SUBL 3

#include <string>

class lattice
{

public:
  lattice(int l);
  ~lattice();

  int j(int* n, int q);  //get the linear lattice position
  
  void getnq(int*, int&, int);
  void torus(int*); //bringing back to torus coordinates using boundary conditions
  //void torus();

//  int nx, ny, q; //coordinates uniquely labelling a point in the lattice
  
//  int* neighbors;
//  void findn (int, int);
//  void findn0(int, int);

  int*** connectivity; //connectivity matrix of the lattice (all links from a given site)
  int*** links; //half of the links from a given site such that no double counting occurs

  int*** plaquette3; //three-site plaquettes
  int*** plaquette4; //four-site plaquettes

  void print();

//some functions to generate the connectivity matrices for the lattices
  void set_square();
  void set_triangular();
  void set_checkerboard();
  void set_kagome();

  void set_chain();

  void adjacency(int);
  std::string get_desc();

private:
  int L, Q;
  int LD, N;

  int nbr; //number of neighbors stored in the connectivity matrix

  std::string desc;
  int *n1, *n2;
};

#endif
