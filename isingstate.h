

//The class isingstate contains a general ising state on a lattice

#ifndef ISINGSTATE_H
#define ISINGSTATE_H

using namespace std;

/*#include <stdlib.h>
#ifndef GET_RAND
#define GET_RAND ( (double)rand()/((double)RAND_MAX + 1.0) ) // the standard c library in stdlib.h
#endif

#include "randoma.h"
*/

#include "helperfunctions.h"
#include "lattice.h"

//number of states per lattice site
#define NS 2

class isingstate
{

public:
  //Linear system size L for LxL square lattice
  isingstate(int L); //trivial unit cell
  isingstate(int L, int ucell); // non-trivial unit cell
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

  //flavor number per site
//  int na;

  int N; //N is the number of sites

  int* Na;

  lattice* mylattice;

  void init(int);
  void init(int, int);// version for the case of non-trivial unit cell size
  int save(const std::string);
  int load(const std::string);
};


#endif

