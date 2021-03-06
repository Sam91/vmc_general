#include "isingstate.h"
#include "mysql_wrapper.h"
#include <iomanip>

//This is a virtual class for a wave function

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#define WFC 0
#define JASTROW 0

//number of real observables
//#define NO 5 //this is now set in wavefunction.cpp (bcs we want to make it depend on system size)

//number of variational parameters for which we would like to calculate the gradient (not implemented)
#define NP 0

#define SMALL 1e-14

struct baseparameters
{
  int N; //Number of sites/factors in the local tensor product

  string desc; //an english description of the state

  bool *ap; //spinon boundary conditions
};

class wavefunction
{

public:
  wavefunction(int);
  virtual ~wavefunction();

  virtual void getwf() =0; //calculate <alpha|wf> from scratch
#if WFC
  virtual complex<double> swap(int i1, int i2, bool ratio) =0; //swap states on two sites and return the wf ratio (for ratio=true, no backup is needed)
  virtual complex<double> crop(int i1, int i2, bool ratio) =0; //particle non-conserving process and return the wf ratio
#else
  virtual double swap(int i1, int i2, bool) =0; //swap states on two sites and return the wf ratio
  virtual double crop(int i1, int i2, bool) =0; //particle non-conserving process and return the wf ratio
#endif

  virtual bool step(); //do a MC step in the ising space and accept it with probability (w'/w)^2; overwrite for flavor non-conserving moves.

  virtual void find_starting_conf();

  virtual void backup_data() =0;
  virtual void restore_data() =0;

  virtual void print() =0;
  virtual int insert_db() =0;
  virtual void print_avgs();
  //virtual int insert_file(const char*) =0;
  virtual int insert_file(string) =0;

  int getL();
  int getN();
  void print_alpha();
  void set_random_conf();

  void walk(); //a random walk that should generate an independent state in sim L2 steps
  void set_walk_length(int);
  int get_walk_length();

  void set_mc_length(int);  //set the length of MC over which we accumulate date

  void accumulate();
  void walk_accumulate();
  void collect_data();
  void calculate_statistics();

  void initiate_f( int );
  void destroy_f( int );

  //this routine calculates the the wf gradient for the parameters specified
  //virtual void get_dwf() =0;

#if WFC
  complex<double> wf, wf_old;
  complex<double> *dwf;
#else
  double wf, wf_old;
  double *dwf;
#endif

  int accepted;

  baseparameters* bpars;
  
  double *js; //jastrow factors

  double jl; //dogleg term (triangular lattice)

  void set_lattice(const string&);

  int save_alpha();
  int load_alpha();

  virtual void correct_cff(bool) =0; //corrects the wave function normalization by a factor
  void reset_run();

  void calculate_exact();
  void accumulate_exact();
  void print_f0();
  void print_NF();

protected:
  isingstate* alpha;

  int L, Q, LD, N; //system size

  int NO; //number of operators to measure

  bool accept( complex<double> );
  bool accept( double );

  int run;
  int walkl; //walk length
  int mc_length;
  int nk; //number of bins

  double jastrow(int,int,double,double);
  virtual double jastrow();
  virtual double jastrow(int,int);
  //special jastrows for the anisotropic triangular lattice
  virtual double jastrow3();
  virtual double jastrow3(int,int);

  virtual complex<double> dogleg(); //dogleg functions for the triangular lattice
  virtual complex<double> dogleg(int, int);

  double *average, *sigma;
  double **daverage, **dsigma; //gradient terms

  double tPi, tPiL;
  complex<double> I;

  int* NF; //Flavor number to which this wf is projected
  double cff;
  double norm;

private:
  double **f, ***fd;
  double *f0, **f0d;
  double *fj;
};

#endif
