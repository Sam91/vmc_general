#include "paired2k.h"

#ifndef AMPEREAN_H
#define AMPEREAN_H

//a structure containing relevant parameters of the wave function
struct parameters : public baseparameters
{
  double t1; //nearest-neighbor hopping
  double t2; //next-neighbor

  double mu; //chemical potential (same for both flavors)

  int* NF; //this wave function is a flavor number eigenstate, and these are the flavor numbers

  double dd,dd0; //maximal paring gap
  double lth, lr1, lr2;

  int r; //
};

class amperean : public paired2k
{
public:
  amperean(int, int);
  ~amperean();

  virtual void create_dd(); //sets delta, xi, mu, ap in paired2k
  virtual int insert_db();
  virtual void print();

  void findmu();
  void findmu(int);

  parameters* pars;

  void get_polar(double&, double&, int, int);
  void get_carth(int&, int&, double, double);

  void get_fflo(int &k1, int &k2, int j1, int j2, int r);
  static double xif(double, double, double);
  int getQ(int, int, int);
  double getQ_angle(int, int);
private:
  //double theta;
  double eps;
  double xik(double, double);
};

#endif
