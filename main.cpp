
#include <iostream>

using namespace std; 

#include "isingstate.h"

int fact(int n)
{
  int r = 1;
  for(int i=1; i<=n; i++)
    r *= i;
  return r;
}

//int main()
int main(int argc, char *argv[])
{
  //cout << "Hello world\n";

  int na = 3;
  int L = 12;

  isingstate alpha(na, L);

  int *Na = new int[na];
  Na[0] = 48; Na[1] = 48;

  alpha.set_random_conf( Na );
  //alpha.print();

  int *hist = new int[L*L/2+1];

  for(int k=0; k<120000; k++) {
    alpha.set_random_conf( Na );
    //cout << "Na[1]: " << alpha.getNa(1) << "\n";

    for(int i=0; i<400; i++)
      alpha.step();
    //cout << "Na[1]: " << alpha.getNa(1) << "\n";

    hist[alpha.getNa(1)]++;
  }

  alpha.print();

  cout << "Na test: " << alpha.checkNa() << "\n";

  for(int k=0; k<=L*L/2; k++)
    cout << "(" << k << "," << hist[k] << "),\n";
  cout << "\n";

//  for(int k=0; k<L*L/2; k++)
//    cout << fact(L*L)/(fact(k)*fact(k)*fact(L*L-2*k)) << ", ";
//  cout << "\n";

  save_v(hist, L*L/2+1, "hist.out");

  return 1;
}

