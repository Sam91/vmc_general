#include "wavefunction.h"
#include <iomanip>

wavefunction::wavefunction(int n, int l) //arguments: number of sites and linear system size
{
  this->N = n;
  cout << "WF: setting sites to " << n << "\n";

  this->L = l; this->L2 = l*l; //number of unit cells (number of sites depends on the lattice)

  alpha = new isingstate( L );

  f0 = new double[NO];
  for(int no=0; no<NO; no++) f0[no] = 0.;

  fj = new double[NO];

#if NP
  //gradient vector; for each operator and for the wave function itself
  daverage = createdouble(NO, NP);
  dsigma = createdouble(NO, NP);

  f0d = createdouble(NO+1,NP);
  for(int no=0; no<=NO; no++)
    for(int np=0; np<NP; np++) f0d[no][np] = 0.;

#if WFC
  dwf = new complex<double>[NP];
#else
  dwf = new double[NP];
#endif
#endif

  average = new double[NO];
  sigma   = new double[NO];
  
  NF = new int[NS];
  for(int n=0; n<NS; n++) NF[n] = N/NS;

  tPiL = 2.*M_PI/(double)L;
  tPi  = 2.*M_PI;
  I = complex<double>(0.,1.);
  cff = 1.;
  
  //jastrow factors
  js = new double[10];
  for(int i=0; i<10; i++) js[i] = 0.;

  //dogleg variational parameter
  jl = 0.;

  //default mc length
  mc_length = 40;
}

wavefunction::~wavefunction()
{
  cout << "wavefunction::~wavefunction()\n";
  delete[] f0;
  delete[] fj;

  delete[] average;
  delete[] sigma;

  delete[] NF;

  delete alpha;
#if NP
  destroy(daverage, NO);
  destroy(dsigma, NO);
  delete[] dwf;
  destroy(f0d, NO);
#endif
  delete[] js;
}

void wavefunction::reset_run() {run = 0;}

void wavefunction::set_lattice(const string& s)
{
  if( s.compare("square")==0 )
    alpha->mylattice->set_square();

  if( s.compare("checkerboard")==0 )
    alpha->mylattice->set_checkerboard();

  if( s.compare("triangular")==0 )
    alpha->mylattice->set_triangular();

  //alpha->mylattice->adjacency( 0 );

  if( bpars->N != alpha->N ) {
    cout << "ERROR: number of sites not consistent!\n";
    exit(-1);
  }
}

void wavefunction::set_random_conf()
{
  int *N0 = new int[NS];
  for(int n=0; n<NS; n++) N0[n] = L2/NS;
  alpha->set_random_conf( N0 );
  delete[] N0;
}

int wavefunction::save_alpha()
{
  return alpha->save("a.dat");
}

int wavefunction::load_alpha()
{
  return alpha->load("a.dat");
}

void wavefunction::print_alpha() { alpha->print(); }

bool wavefunction::accept(double p)
{
  double p2 = p*p;

  if( p2>= 1. ) return true;

  if( GET_RAND < p2 ) return true;
  else return false;
}

bool wavefunction::accept(complex<double> p)
{
  double p2 = pow(abs(p),2);

  //cout << " p2= " << p2 << "\n";

  if( p2 >= 1. ) { /*cout << "  true\n";*/ return true;}

  //double rn = GET_RAND;
  //cout << "rand = " << rn << "\n";

  if( GET_RAND < p2 ) {/*cout << "  true\n";*/ return true;}
  else {/*cout << "  false\n";*/ return false;}
}

int wavefunction::getL() {return L;}

//here we accumulate the observables to be measured
void wavefunction::accumulate()
{
  //cout << "accumulate << " << k << "\n";

  //get the derivative of the wave function and accumulate them
#if NP
  get_dwf();
  for(int p=0; p<NP; p++) f0d[NO][p] += dwf[p]; //derivative of the norm 
#endif

  for(int i=0; i<NO; i++) fj[i] = 0.;

#if (NS==2) && (NO>2) //in this case, we calculate ring terms

  int i2,i3,i4;
#if WFC
  complex<double> r;
#else
  double r;
#endif

#endif

  //double p1, p2, p3;

  if( abs(wf)<SMALL ) {
    cout << "Vanishing wf in accumulate...\n";
    throw 3;
  }

  //accumulate the nearest neighbor swap operator
  for(int i1=0; i1<N; i1++ )
  {
    //alpha->mylattice->findn0(i1, 1);
    //cout << "i1: " << i1 << ", ";
    //for(int i=0; i<5; i++) cout << alpha->mylattice->neighbors[i] << "; ";
    //cout << "\n";

    //nearest neighbors
/*
     for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)
    {
#if WFC
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true ).real();  //swap with ratio=true(i.e., only the ratio is calculated but the wf is not updated)
#else
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true );
#endif
    }
*/
    //works for triangular lattice only

#if WFC
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][1], true ).real();
      fj[4] += swap( i1, alpha->mylattice->links[i1][0][2], true ).real();
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][3], true ).real();
#else
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][1], true );
      fj[4] += swap( i1, alpha->mylattice->links[i1][0][2], true );
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][3], true );
/*
      p1 = swap( i1, alpha->mylattice->links[i1][0][1], true );
      p2 = swap( i1, alpha->mylattice->links[i1][0][2], true );
      p3 = swap( i1, alpha->mylattice->links[i1][0][3], true );
      fj[0] += p1;
      fj[4] += p2;
      fj[0] += p3;
      fj[4+3*i1+1] += p1;
      fj[4+3*i1+2] += p2;
      fj[4+3*i1+3] += p3;
*/
#endif


    //next-nearest neighbor exchange
    for(int n=1; n<=alpha->mylattice->links[i1][1][0]; n++)
    {
#if WFC
      fj[1] += swap( i1, alpha->mylattice->links[i1][1][n], true ).real();
#else
      fj[1] += swap( i1, alpha->mylattice->links[i1][1][n], true );
#endif
    }

    //third-neighbor exchange
    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)
    {
#if WFC
      fj[2] += swap( i1, alpha->mylattice->links[i1][2][n], true ).real();
#else
      fj[2] += swap( i1, alpha->mylattice->links[i1][2][n], true );
#endif
    }

#if (NS==2) && (NO>2)

    //four-site plaquette term (this measurement is optimized for TWO flavors, where a four-cycle does not arise)
    //loop over plaquettes
//cout << "Looping " << alpha->mylattice->plaquette4[i1][0][0] << " plaquettes.\n";
    for(int n=1; n<=alpha->mylattice->plaquette4[i1][0][0]; n++)
    {
//cout << "n: " << n << "\n";

      i2 = alpha->mylattice->plaquette4[i1][0][3*n-2];
      i3 = alpha->mylattice->plaquette4[i1][0][3*n-1];
      i4 = alpha->mylattice->plaquette4[i1][0][3*n-0];

      //cout << "Plaquette: " << "(" << i1 << "," << i2 << "," << i3 << "," << i4 << ")\n";

      if( (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(aaaa)
      {
        fj[3] += 1.;
        continue;
      }

      if( (alpha->lconf[i1]!=alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(abbb)
      {
//        r = swap(i1, i2, true);
//        if(abs(r)>50.) cout << "ERROR r1a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i1, i2, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i1, i2, true);
#endif
        continue;
      }

      if( (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) ) //(babb)
      {
//        r = swap(i2, i3, true);
//        if(abs(r)>50.) cout << "ERROR r1b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i2, i3, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i2, i3, true);
#endif
        continue;
      }

      if( (alpha->lconf[i3]!=alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) && (alpha->lconf[i1]==alpha->lconf[i2]) ) //(bbab)
      {
//        r = swap(i3, i4, true);
//        if(abs(r)>50.) cout << "ERROR r1c, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i3, i4, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i3, i4, true);
#endif
        continue;
      }

      if( (alpha->lconf[i4]!=alpha->lconf[i1]) && (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) ) //(bbba)
      {
//        r = swap(i4, i1, true);
//        if(abs(r)>50.) cout << "ERROR r1d, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i4, i1, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i4, i1, true);
#endif
        continue;
      }

      if( (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(aabb)
      {
//        r = swap(i1, i3, true);
//        if(abs(r)>50.) cout << "ERROR r2a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i1, i3, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i1, i3, true);
#endif
        continue;
      }

      if( (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]!=alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) ) //(abba)
      {
//        r = swap(i2, i4, true);
//        if(abs(r)>50.) cout << "ERROR r2b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i2, i4, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i2, i4, true);
#endif
        continue;
      }

      if( (alpha->lconf[i1]!=alpha->lconf[i2]) && (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]!=alpha->lconf[i4]) ) //(abab); this is a product of two transpositions
      {
        backup_data();
        int c;

        r = swap(i1, i2, false); //here, we update the matrix
//        if(abs(r)>50.) cout << "ERROR r4a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";

        if( abs(r)<SMALL ) //intermediate wf vanishes
        {
          c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
          r = wf_old;
          getwf();
//          r = wf/r;
//          if(abs(r)>50.) cout << "ERROR r4a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//          fj[2] += r.real();
          fj[3] += (wf/r).real();
#else
//          fj[2] += r;
          fj[3] += wf/r;
#endif
          restore_data();
          c = alpha->lconf[i1]; alpha->lconf[i1]=alpha->lconf[i2]; alpha->lconf[i2]=c;
          c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
          continue;
        }

        r *= swap(i3, i4, false); //here, we update the matrix
//        if(abs(r)>50.) cout << "ERROR r4b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";

#if WFC
        fj[3] += r.real();
#else          
        fj[3] += r;
#endif    
        restore_data();
        c = alpha->lconf[i1]; alpha->lconf[i1]=alpha->lconf[i2]; alpha->lconf[i2]=c;
        c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
      }
    } //plaquette loop
#endif //NS==2

  } //end of loop over sites

  for(int no=0; no<NO; no++) {//operator index
    f0[no] += fj[no];

#if NP //gradient calculation for each operator
    for(int p=0; p<NP; p++) f0d[no][p] += fj[no]*dwf[p]; //parameter index
#endif
  }
}

/*void wavefunction::setf(void* ff)
{
  f = (double**)ff;
}*/

void wavefunction::initiate_f( int n )
{
  f = createdouble(n, NO);
#if NP
  fd = new (double(**[n]));
  for(int i=0; i<n; i++) fd[i] = createdouble(NO, NP);
#endif
  run = 0;
  nk = n;
}
void wavefunction::destroy_f( int n )
{
  destroy(f, n);
#if NP
  for(int i=0; i<n; i++) destroy(fd[i], NO);
  delete[] fd;
#endif
}

void wavefunction::collect_data()
{
  for(int no=0; no<NO; no++)
  {
    f[run][no] = f0[no]/((double)(L2*mc_length));
    f0[no] = 0.;
#if NP
    for(int p=0; p<NP; p++) {
      fd[run][no][p] = f0d[no][p]/((double)(L2*mc_length)) - f0d[NO][p]/((double)(mc_length));
      f0d[no][p] = 0.;
    }
    for(int p=0; p<NO; p++) f0d[NO][p] = 0.;
#endif
  }

  run++;
}

void wavefunction::set_walk_length(int wl) { walkl = wl; }

void wavefunction::set_mc_length(int mc) {cout << "Setting mc_length to " << mc << "\n"; mc_length = mc;}

int wavefunction::get_walk_length() {return walkl;}

void wavefunction::walk()
{
  for(int i=0; i<walkl; i++) step();
}

void wavefunction::walk_accumulate()
{
  //cout << "walkl = " << walkl << "\n";
  for(int j=0; j<mc_length; j++) {
    walk();
    accumulate();
  }
}

void wavefunction::calculate_statistics()
{
  for(int no=0; no<NO; no++) {
    average[no] = 0.; sigma[no] = 0.;
#if NP
    for(int p=0; p<NP; p++) { daverage[no][p] = 0.; dsigma[no][p] = 0.; }
#endif
  }
  cout << "runs: " << run << "\n";

  for(int no=0; no<NO; no++)
  {
    for(int i=0; i<run; i++)
    {
      //cout << "f[" << i << "][" << j << "] = " << f[i][j] << "\n";
      average[no] += f[i][no]; sigma[no] += pow(abs(f[i][no]),2);
#if NP
      for(int p=0; p<NP; p++) { daverage[no][p] += fd[i][no][p]; dsigma[no][p] = pow(fd[i][no][p],2); }
#endif
    }
    if(run>1) {
      average[no] /= (double)run; sigma[no] /= (double)run; sigma[no] -= pow(abs(average[no]),2);;
      if( sigma[no]>=0. ) sigma[no] = sqrt( abs(sigma[no]/(run-1)) );
      else {
        sigma[no] = -1.; cout << "Variance issue\n";
      }
#if NP
      for(int p=0; p<NP; p++) {
        daverage[no][p] /= (double)run; dsigma[no][p] /= (double)run; dsigma[no][p] -= pow(abs(daverage[no][p]),2);;
        if( dsigma[no][p]>=0. ) dsigma[no][p] = sqrt( abs(dsigma[no][p]/(run-1)) );
        else {
          dsigma[no][p] = -1.; cout << "Variance issue\n";
        }
      }
#endif
    }
  }
  run = 0;

  for(int no=0; no<NO; no++) {
    cout << std::fixed << setprecision(4) << average[no] << " pm " << sigma[no] << "\n";
#if NP
    for(int p=0; p<NP; p++)
      cout << std::fixed << setprecision(4) << daverage[no][p] << " pm " << dsigma[no][p] << "\n";
#endif
  }
  cout << "\n";
}

/*
 * This basic implementation of the MC step function performs a simple exchange of flavors (swap).
 * It should be overwritten when flavor non-conserving moves should be performed (crop).
*/
bool wavefunction::step()
{
  int i1, i2;
#if WFC
  complex<double> p;
#else
  double p;
#endif

  do {

    //pick two sites at random
    double fr = N*GET_RAND;

    i1 = (int)(fr);
    i2 = (int)( (fr-i1)*N );

    if( alpha->lconf[i1]==alpha->lconf[i2] ) continue; //continue until we find two sites with different states

    backup_data();

    p = swap(i1, i2, false); //this may execute getwf() and updates wf in any case (for true); it also swaps lconf

    if( abs(p)>1e200 )
    {
      cout << "ERROR: divergent wf\n";
      exit(-1);
    }

//    getwf(); complex<double> pp = wf/wf_old;
//    if( (abs(pp-p)>1e-2) ) {
//      cout << "Warn: p = " << p << "; pp = " << pp << ". i1 = " << i1 << "; i2 = " << i2 << "\n";
//    }
//    else cout << "good agreement\n";

    if( accept(p) )
    {
      accepted++;
      //cout << "accepted: " << std::scientific << setprecision(1) << wf << "; |p|2=" << std::fixed << abs(p)*abs(p) << "\n";
  
      return true;
    } else
    {
      //cout << "rejected: " << std::scientific << setprecision(1) << wf << "; |p|2=" << std::fixed << abs(p)*abs(p) << ")\n";
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c; //swap back the states in alpha
      restore_data();

      return false;
    }
  } while (true);
}

void wavefunction::find_starting_conf()
{
  int* Na0 = new int[NS];

  for(int n=0; n<NS; n++) Na0[n] = NF[n];

  //find_max_conf();
/*  int i, imax = 10000;
  for(i=0; i<imax; i++)
  {
    alpha->set_random_conf( Na0 );
    getwf();
    if( abs(wf)>SMALL ) break;
  }
  if( i==imax) cout << "ERROR: cannot find non-zero configuration...\n";
*/

  alpha->set_random_conf( Na0 );
  getwf();

  //cout << "Max config:\n";
  //alpha->print();
  //cout << "wf: " << wf << "\n";

  delete[] Na0;
}

//find the jastrow coefficient between two sites
double wavefunction::jastrow(int i1, int i2, double j1, double j2)
{
  double r;
  if( (alpha->lconf[i1]==1 &&  alpha->lconf[i1]==1) || (alpha->lconf[i1]==2 &&  alpha->lconf[i1]==2) ) r = j1;
  else if( (alpha->lconf[i1]==1 &&  alpha->lconf[i1]==2) || (alpha->lconf[i1]==2 &&  alpha->lconf[i1]==1) ) r = -j1;
  else r=1.;

  if( alpha->lconf[i1]!=0 && alpha->lconf[i2]!=0 ) return r+j2;
  else return r;
}

//accumulate the nearest neighbor Jastrow factor for this state
double wavefunction::jastrow()
{
  int nJ1=0;
  int nJ2=0;
  int nJ3=0;
  
  for(int i1=0; i1<N; i1++)
  {
    //accumulate permutation invariant Jastrow factors
    for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)//nearest-neighbor
    {
      if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][n]] ) nJ1 += 2; //add two when same flavor occupies
      nJ1--; //subtract 1 per link
    }
    
    for(int n=1; n<=alpha->mylattice->links[i1][1][0]; n++)//next-neighbor
    {
      if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][n]] ) nJ2 += 2;
      nJ2--;
    }

    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)//third-neighbor
    {
      if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][2][n]] ) nJ3 += 2;
      nJ3--;
    }
  }
  
  double r = exp( - (js[0]*nJ1 + js[1]*nJ2 + js[2]*nJ3)/2. );
  if( abs(r)>1e290 ) {
    cout << "diverging jastrow\n";
    throw 20;
  }
  cout << "JF0 = " << r << "\n";
  return r;
}

double wavefunction::jastrow(int i1, int i2 )
{
  if( alpha->lconf[i1] == alpha->lconf[i2] ) return 1.;

  int nJ1=0;
  int nJ2=0;
  int nJ3=0;

  //loop over the nearest neighbors of i1
  for(int n=1; n<=alpha->mylattice->connectivity[i1][0][0]; n++)
  {
    if( i2 == alpha->mylattice->connectivity[i1][0][n] ) continue;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][n]] ) nJ1--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][n]] ) nJ1++;
  }

  //loop over the nearest neighbors of i2
  for(int n=1; n<=alpha->mylattice->connectivity[i2][0][0]; n++)
  {
    if( i1 == alpha->mylattice->connectivity[i2][0][n] ) continue;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][n]] ) nJ1--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][n]] ) nJ1++;
  }

  //loop over the next neighbors of i1
  for(int n=1; n<=alpha->mylattice->connectivity[i1][1][0]; n++)
  {
    if( i2 == alpha->mylattice->connectivity[i1][1][n] ) continue;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][n]] ) nJ2--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][n]] ) nJ2++;
  }

  //loop over the next neighbors of i2
  for(int n=1; n<=alpha->mylattice->connectivity[i2][1][0]; n++)
  {
    if( i1 == alpha->mylattice->connectivity[i2][1][n] ) continue;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][n]] ) nJ2--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][n]] ) nJ2++;
  }

  //loop over the third neighbors of i1
  for(int n=1; n<=alpha->mylattice->connectivity[i1][2][0]; n++)
  {
    if( i2 == alpha->mylattice->connectivity[i1][2][n] ) continue;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3++;
  }

  //loop over the third neighbors of i2
  for(int n=1; n<=alpha->mylattice->connectivity[i2][2][0]; n++)
  {
    if( i1 == alpha->mylattice->connectivity[i2][2][n] ) continue;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3++;
  }

  double r = exp( -(js[0]*nJ1 + js[1]*nJ2 + js[2]*nJ3) );
  
  //cout << "JF = " << r << "\n";
  return r;
}


//jastrow factors for the anisotropic triangular lattice

//accumulate the nearest neighbor Jastrow factor for this state

double wavefunction::jastrow3()
{
  int nJ1a=0;
  int nJ1b=0;
  int nJ1c=0;
  int nJ2a=0;
  int nJ2b=0;
  int nJ2c=0;
  int nJ3=0;
  
  for(int i1=0; i1<N; i1++)
  {
    //accumulate permutation invariant Jastrow factors
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][1]] ) nJ1a += 2; //add two when same flavor occupies
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][2]] ) nJ1b += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][0][3]] ) nJ1c += 2;

    nJ1a-=1; //subtract 1 per link
    nJ1b-=1;
    nJ1c-=1;

    //next neighbor
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][2]] ) nJ2a += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][1]] ) nJ2b += 2;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][1][3]] ) nJ2c += 2;

    nJ2a-=1;
    nJ2b-=1;
    nJ2c-=1;

    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)//third-neighbor
    {
      if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->links[i1][2][n]] ) nJ3 += 2;
      nJ3--;
    }
  }
  
  double r = exp( - (js[0]*nJ1a + js[1]*nJ2a + js[2]*nJ3 + js[3]*nJ1b + js[4]*nJ2b + js[5]*nJ1c + js[6]*nJ2c)/2. );
  if( abs(r)>1e290 )
  {
    cout << "diverging jastrow\n";
    throw 20;
  }
  for(int i=0; i<7; i++) cout << js[i] << ", "; cout << "\n";
  cout << "JF0 = " << r << "\n";
  return r;
}

double wavefunction::jastrow3(int i1, int i2 )
{
  if( alpha->lconf[i1] == alpha->lconf[i2] ) return 1.;

  int nJ1a=0;
  int nJ1b=0;
  int nJ1c=0;
  int nJ2a=0;
  int nJ2b=0;
  int nJ2c=0;
  int nJ3=0;

  //loop over the nearest neighbors of i1
  if( i2 != alpha->mylattice->connectivity[i1][0][1] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][1]] ) nJ1a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][1]] ) nJ1a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][2] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][2]] ) nJ1b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][2]] ) nJ1b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][3] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][3]] ) nJ1c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][3]] ) nJ1c++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][4] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][4]] ) nJ1a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][4]] ) nJ1a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][5] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][5]] ) nJ1b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][5]] ) nJ1b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][0][6] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][0][6]] ) nJ1c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][0][6]] ) nJ1c++;
  }

  //loop over the nearest neighbors of i2
  if( i1 != alpha->mylattice->connectivity[i2][0][1] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][1]] ) nJ1a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][1]] ) nJ1a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][2] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][2]] ) nJ1b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][2]] ) nJ1b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][3] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][3]] ) nJ1c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][3]] ) nJ1c++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][4] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][4]] ) nJ1a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][4]] ) nJ1a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][5] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][5]] ) nJ1b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][5]] ) nJ1b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][0][6] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][0][6]] ) nJ1c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][0][6]] ) nJ1c++;
  }

  //loop over the next neighbors of i1
  if( i2 != alpha->mylattice->connectivity[i1][1][1] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][1]] ) nJ2a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][1]] ) nJ2a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][2] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][2]] ) nJ2b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][2]] ) nJ2b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][3] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][3]] ) nJ2c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][3]] ) nJ2c++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][4] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][4]] ) nJ2a--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][4]] ) nJ2a++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][5] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][5]] ) nJ2b--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][5]] ) nJ2b++;
  }
  if( i2 != alpha->mylattice->connectivity[i1][1][6] )
  {
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][1][6]] ) nJ2c--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][1][6]] ) nJ2c++;
  }

  //loop over the next neighbors of i2
  if( i1 != alpha->mylattice->connectivity[i2][1][1] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][1]] ) nJ2a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][1]] ) nJ2a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][2] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][2]] ) nJ2b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][2]] ) nJ2b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][3] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][3]] ) nJ2c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][3]] ) nJ2c++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][4] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][4]] ) nJ2a--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][4]] ) nJ2a++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][5] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][5]] ) nJ2b--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][5]] ) nJ2b++;
  }
  if( i1 != alpha->mylattice->connectivity[i2][1][6] )
  {
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][1][6]] ) nJ2c--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][1][6]] ) nJ2c++;
  }

  //loop over third neighbors of i1
  for(int n=1; n<=alpha->mylattice->connectivity[i1][2][0]; n++)
  {
    if( i2 == alpha->mylattice->connectivity[i1][2][n] ) continue;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3--;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i1][2][n]] ) nJ3++;
  }

  //loop over third neighbors of i2
  for(int n=1; n<=alpha->mylattice->connectivity[i2][2][0]; n++)
  {
    if( i1 == alpha->mylattice->connectivity[i2][2][n] ) continue;
    if( alpha->lconf[i2] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3--;
    if( alpha->lconf[i1] == alpha->lconf[alpha->mylattice->connectivity[i2][2][n]] ) nJ3++;
  }

  double r = exp( -(js[0]*nJ1a + js[1]*nJ2a + js[2]*nJ3 + js[3]*nJ1b + js[4]*nJ2b + js[5]*nJ1c + js[6]*nJ2c) );
  
  //cout << "JF = " << r << "\n";
  return r;
}

// Calculates the dogleg factor for the triangular lattice
complex<double> wavefunction::dogleg()
{
  int l = 0;
  
  for(int i1=0; i1<N; i1++) //loop through the middle point of the dog leg
  {
//    if( ((i1%L)+(i1/L))%2 == 0 )
//    {
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ]-2) == 1) l++; else l--;
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ]-2) == 1) l--; else l++;
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ]-2) == 1) l++; else l--;
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ]-2) == 1) l--; else l++;
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ]-2) == 1) l++; else l--;
      if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ]-2) == 1) l--; else l++;
/*    } else
    {
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ])-3) == 3) l--; else l++;
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ])-3) == 3) l++; else l--;
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ])-3) == 3) l--; else l++;
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ])-3) == 3) l++; else l--;
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ])-3) == 3) l--; else l++;
      if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ])-3) == 3) l++; else l--;
    }
*/
  }
//  cout << "dogleg(): l = " << l << endl;
  return exp( complex<double>(0.,(double)l * (this->jl) ) );
}

//returns the change in dogleg after as swap
complex<double> wavefunction::dogleg(int i1, int i2)
{
  int l = 0;

  if( alpha->lconf[i1]==alpha->lconf[i2] ) return 1.;

//  if( ((i1%L)+(i1/L))%2 == 0 ) //dogleg for i1
//  {
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ]-2) == 1) l++; else l--; //add for the new state
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ]-2) == 1) l--; else l++;

    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ]-2) == 1) l--; else l++; //remove the old one
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ]-2) == 1) l++; else l--;

/*  } else
  {
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ])-3) == 3) l++; else l--;

    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][1] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i1][0][2] ])-3) == 3) l--; else l++;
  }
*/

//  if( ((i2%L)+(i2/L))%2 == 0 ) //dogleg for i2
//  {
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ]-2) == 1) l++; else l--; //add for the new state
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ]-2) == 1) l--; else l++;

    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ]-2) == 1) l--; else l++; //remove the old one
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ]-2) == 1) l++; else l--;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ]-2) == 1) l--; else l++;
    if( abs(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ]-2) == 1) l++; else l--;

/*  } else
  {
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i1] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ])-3) == 3) l++; else l--;

    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][3] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][4] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ])-3) == 3) l--; else l++;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][5] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][1] ])-3) == 3) l++; else l--;
    if( abs(2*(alpha->lconf[i2] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][6] ] + alpha->lconf[ alpha->mylattice->connectivity[i2][0][2] ])-3) == 3) l--; else l++;
  }
*/
//  cout << "dogleg(i1,i2): l = " << l << endl;
  return exp( complex<double>(0., (double)l * (this->jl) ) );
}
