//The class isingstate contains a general ising state on a lattice

#ifndef ISINGSTATE_H
#define ISINGSTATE_H

using namespace std;


#include "helperfunctions.h"
#include "lattice.h"

//number of states per lattice site
#define NS 2

class isingstate
{

public:
  //Linear system size L for LxL square lattice
  isingstate(int L); //trivial unit cell
  //isingstate(int L, int ucell); // non-trivial unit cell
  isingstate();

  ~isingstate();

  //set a random conf with total flavor numnber N
  void set_random_conf(int *N);

  //print the current configuration
  void print();

  bool checkNa();
  int getNa(int);

  void step();

  int* lconf;

  int N; //N is the number of sites

  int* Na;

  lattice* mylattice;

  void init(int);
  //void init(int, int);// version for the case of non-trivial unit cell size
  int save(const std::string);
  int load(const std::string);

  void iter(); //iterate to the next ising state
  int getint(); //get the integer representation
  void setfirst(); //set the lowest state (int rep)
  bool islast(); //check if it is the last state
};


#endif

