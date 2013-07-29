
#include "randoma.h"
#include "isingstate.h"

//na : number of states per site
//sys: linear size of the 2D bravais lattice
//q  : number of sites per sublattice
isingstate::isingstate(int sys)
{
  init(sys);
}

isingstate::isingstate()
{
  init(4);
}

isingstate::~isingstate()
{
  cout << "Destroying Ising state\n";
  destroy(conf, L);
  delete[] Na;
  delete[] lconf;
  delete mylattice;
}

void isingstate::init( int sys )
{
  cout << "Creating new Ising state with NS=" << NS << " and L=" << sys << "\n";

  //this->na = na; 
  this->L = sys;
  this->N = L*L;

  conf = createint(L);
  lconf = new int[N];
  Na = new int[NS];

  mylattice = new lattice(sys, 0);

  randomize();
}

int isingstate::save(const std::string s) //save the configuration to a file
{
  return save_v(lconf, N, s.c_str() );
}

int isingstate::load(const std::string s)
{
  return load_v(lconf, N, s.c_str() );
}

void isingstate::print()
{
  cout << "Na = ";
  for(int n=0; n<NS; n++) cout << Na[n] << ", ";
  cout << "\n";
  cout << "conf:\n";
  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++) {
      cout << lconf[i+L*j] << ", ";
    }
    cout << "\n";
  }
  cout << "\n";
}

bool isingstate::checkNa()
{
  int* Natest = new int[NS];
  for(int i=0; i<NS; i++) Natest[i] = 0;

  for(int i=0; i<N; i++)
    //for(int j=0; j<L; j++)
      Natest[ lconf[i] ]++;

  for(int i=0; i<NS; i++)
    if( Natest[i]!=Na[i] ) {
      delete[] Natest;
      return false;
    }

  delete[] Natest;
  return true;  
}

int isingstate::getNa(int k){ return Na[k]; }

void isingstate::set_random_conf(int *Na0)
{
  int n = 0;
  for(int i=0; i<NS-1; i++) {
    if(Na0[i]<0) {
      cout << "Error in flavor number "<< i << "\n";
      return;
    }
    n += Na0[i];
  }
  if( n>N ) {
    cout << "Error in total flavor number\n";
    return;
  }
  Na0[NS-1] = N-n;

  int *ourNa = new int[NS];
  for(int i=0; i<NS; i++) ourNa[i] = Na[i] = Na0[i];
  
  int N0 = N;
  //double *Nap = new double[na];
  int Nr, Nat;

  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++) {
      //for(int k=0; k<na; k++)
      //  Nap[k] = (double)Na0[k]/(double)N0;

      Nr = (int)((double)N0*GET_RAND);
      Nat = 0;

      //cout << "N0=" << N0 << "\n";
      //cout << "Nr=" << Nr << "\n";

      for(int k=0; k<NS; k++) { //conf = 0 ... NS-1
        Nat += ourNa[k];

        if( Nr<Nat ) { //found a configuration
          conf[i][j] = k;
          lconf[i+L*j] = k;
          ourNa[k]--; N0--;
          break;
        }
      }
    
    }
  }
  delete[] ourNa;
}

/*
//make an elementary step inside the constraint space
//this now depends on the particular hilbert space we want to consider

//This generates the space na=3 with the constraint N1 = N2 and 00->12,21; 12->00,21; 21->00,12
void isingstate::step()
{
  int i1, i2;
  do {
 
    //choose two sites at random
    double fr = N*GET_RAND;

    i1 = (int)(fr);
    i2 = (int)((fr-i1)*N);

    if( i1==i2 ) continue;

    if( lconf[i1]==lconf[i2] )
    {
      if( lconf[i1]==0 ) {
//        fr = GET_RAND;
//        if(fr<.5) { lconf[i1] = 1; lconf[i2] = 2; Na[0]-=2; Na[1]++; Na[2]++;}
//        else { lconf[i1] = 2; lconf[i2] = 1; Na[0]-=2; Na[1]++; Na[2]++; }

        //fr = GET_RAND;
        //if( fr > ((double)((Na[1]+1)*(Na[1]+1))/(double)(Na[0]*(Na[0]-1))) ) continue;
        //if( fr>.5 ) break;
        //if( fr > 1./((double)((Na[0]-1)*(Na[0]-1))/(double)(N*N)) ) continue;

        lconf[i1] = 1; lconf[i2] = 2; Na[0]-=2; Na[1]++; Na[2]++;
        break;
      } else //continue; //the same state; take different sites
             break; //accept the state
    } else {
      if( lconf[i1]==1 && lconf[i2]==2 ) {
        fr = GET_RAND;
        if(fr<.5) { lconf[i1] = 2; lconf[i2] = 1;}
        else { 
          //fr = GET_RAND;
          //if( fr > 1./((double)((Na[0]+1)*(Na[0]+2))/(double)(Na[1]*Na[1])) ) continue;
          lconf[i1] = 0; lconf[i2] = 0; Na[0]+=2; Na[1]--; Na[2]--;
        }
        break;
      }
      if( lconf[i1]==2 && lconf[i2]==1 ) {
        fr = GET_RAND;
        if(fr<.5) { lconf[i1] = 1; lconf[i2] = 2;}
        else {
          //fr = GET_RAND;
          //if( fr > 1./((double)((Na[0]+1)*(Na[0]+2))/(double)(Na[1]*Na[1])) ) continue;
          lconf[i1] = 0; lconf[i2] = 0; Na[0]+=2; Na[1]--; Na[2]--;
        }
        break;
      }
      //just exchange the states
      int c = lconf[i1]; lconf[i1] = lconf[i2]; lconf[i2] = c;
      break;
    }
  } while (true);

  //cout << "Moving (" << i1/L << "," << i1%L << "), (" << i2/L << "," << i2%L << ")\n";
}
*/

