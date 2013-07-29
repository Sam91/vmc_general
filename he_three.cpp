#include "huseelser.h"

//abstract Huse-Else wave function

huseelser::huseelser(int n, int l) : wavefunction( n, l )
{
  //create the d vectors (defining the product state)
#if WFC
  d = createcomplex(NS, N);
#else
  d = createdouble(NS, N);
#endif

}

huseelser::~huseelser()
{
  destroy(d, NS);
}

void huseelser::backup_data()
{
  wf_old = wf;
}

void huseelser::restore_data()
{
  wf = wf_old;
}

//get the wf = <alpha|psi>
void huseelser::getwf()
{
  //product state wave function is given by a complex d = u + iv  vector per site
  //here, we use the 0/up/dn basis
  //u[0...2, site-index]
  //v[0...2, site-index]

  //wf = (double)(L2*L2);
  wf = 1.;

  //cout << "wf norm: " << wf << "\n";

  //loop over all sites and calculate the product state
  for(int i1=0; i1<L2; i1++)
  {
    wf *= 2.*d[alpha->lconf[i1]][i1];
  }

  wf *= jastrow();

  //return wf;
}

//swap operator of two sites i1 and i2 (corresponding to virtual_replacement)
#if WFC
complex<double> huseelser::swap(int i1, int i2)
{
  complex<double> r;
#else
double huseelser::swap(int i1, int i2)
{
  double r;
#endif

  if( alpha->lconf[i1] == alpha->lconf[i2] ) {
    //cout << "identical states in swap\n";
    r = 1.;
  } else {

    if( abs(d[alpha->lconf[i1]][i1])>SMALL && abs(d[alpha->lconf[i2]][i2])>SMALL )
    {
      r = (d[alpha->lconf[i2]][i1]*d[alpha->lconf[i1]][i2])/(d[alpha->lconf[i1]][i1]*d[alpha->lconf[i2]][i2]) * jastrow(i1, alpha->lconf[i2], i2, alpha->lconf[i1]);
    } 
    else {
      //cout << "Hitting node swap\n";
      r = 1.;
    }
  }

  //if( abs(r-1.) > .01 )
  //  cout << "swap (" << alpha->lconf[i1] << ", " << alpha->lconf[i2] << ") returning " << r << "\n";
  return r;

  //additional jastrow terms go here
}

#if WFC
complex<double> huseelser::crop(int i1, int i2) {
  complex<double> r;
#else
double huseelser::crop(int i1, int i2) {
  double r;
#endif

  if( alpha->lconf[i1]==0 && alpha->lconf[i2]==0 ) // 0,0 -> 1,2
  {
    if( abs(d[0][i1])>SMALL && abs(d[0][i2])>SMALL ) {
      r = (d[1][i1]*d[2][i2])/(d[0][i1]*d[0][i2]) * jastrow(i1,alpha->lconf[i2],i2,alpha->lconf[i1]);
    } else {
      //cout << "Hitting node crop\n";
      r = 1.;
    }

  } else { // 1,2 -> 0,0; 2,1 -> 0,0

    if( abs(d[alpha->lconf[i1]][i1])>SMALL && abs(d[alpha->lconf[i1]][i2])>SMALL ) {

      r = (d[0][i1]*d[0][i2])/(d[alpha->lconf[i1]][i1]*d[alpha->lconf[i1]][i2]) * jastrow(i1,alpha->lconf[i2],i2,alpha->lconf[i1]);

    } else {
      cout << "Hitting node crop\n";
      r = 1.;
    }
  }

  //cout << "crop returning " << r << "\n";
  return r;
}

void huseelser::find_starting_conf()
{
  int* Na0 = new int[3];

/*
  for(int n1 = 1; n1<L2/2; n1++) {
    Na0[0] = L2-2*n1; Na0[1] = n1;

    for(int i=0; i<100; i++) {
      alpha->set_random_conf( Na0 );
      cout << "\n";
      alpha->print();
      getwf();
      cout << "wf: " << wf << "\n";
      if( pow(abs(wf),2)>1e-6 ) { delete[] Na0; return; }
    }
  }
  cout << "ERROR: no appropriate starting configuration found\n";
  exit(-1);
*/

  Na0[0] = L2/3;
  Na0[1] = L2/3;

  alpha->set_random_conf( Na0 );

//  find_max_conf();

  getwf();
  //cout << "Max config:\n";
  //alpha->print();
  //cout << "wf: " << wf << "\n";

}

void huseelser::find_max_conf()
{
  double* dabs = new double[NS];
  int nmax; double dmax;

  for(int n=0; n<NS; n++) alpha->Na[n] = 0;

  for(int i=0; i<L2; i++)
  {
    for(int n=0; n<NS; n++)
      dabs[n] = abs(d[n][i]);

    nmax = 0; dmax = dabs[0];
    for(int n=1; n<NS; n++) {
      if( dabs[n]>dmax ) {nmax = n; dmax = dabs[n];}
    }
    alpha->Na[nmax]++;
    //cout << "Nmax = " << nmax << "\n";
    alpha->lconf[i]=nmax;
  }

  delete[] dabs;
}

bool huseelser::step()
{

  int i1, i2;
#if WFC
  complex<double> p;
#else
  double p;
#endif
  do {

    //pick two sites at random
    double fr = (alpha->N)*GET_RAND;

    i1 = (int)(fr);
    i2 = (int)((fr-i1)*alpha->N);

    if( i1==i2 ) continue;

    if( alpha->lconf[i1]==alpha->lconf[i2] ) //identical states on sites
    {
      if( alpha->lconf[i1]==0 ) { //case 0,0 -> 1,-1
/*        fr = GET_RAND;
        if(fr<.5) { alpha->lconf[i1] = 1; alpha->lconf[i2] = 2; alpha->Na[0]-=2; alpha->Na[1]++; alpha->Na[2]++;}
        else { alpha->lconf[i1] = 2; alpha->lconf[i2] = 1; alpha->Na[0]-=2; alpha->Na[1]++; alpha->Na[2]++; }
*/
        //fr = GET_RAND;
        //if( fr > ((double)((alpha->Na[1]+1)*(alpha->Na[1]+1))/(double)(alpha->Na[0]*(alpha->Na[0]-1))) ) continue;
        //if( fr>.5 ) break;
        //if( fr > 1./((double)((alpha->Na[0]-1)*(alpha->Na[0]-1))/(double)(N*N)) ) continue;

        //p = (d[1][i1]*d[2][i2])/(d[0][i1]*d[0][i2]);

        p = crop(i1,i2);
        //possible Jastrows come here

        if( accept(p) ) {
          alpha->lconf[i1] = 1; alpha->lconf[i2] = 2; alpha->Na[0]-=2; alpha->Na[1]++; alpha->Na[2]++;
          backup_data();
          wf *= p;
          accepted++;
          return true;
        }
      } else { //non-0 identical states

        //continue; //the same state; take different sites
        //accepted++; return true; //accept and continue

        do //find two different states and exchange them
        {
          i2 = (int)((alpha->N)*GET_RAND);

        } while( alpha->lconf[i1]==alpha->lconf[i2] || i1==i2);

        p = swap(i1, i2);

        if( accept(p) ) {
          int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
          backup_data();
          wf *= p;
          accepted++; return true;
        }
      }
      break;

    } else { //non-identical states

      if( alpha->lconf[i1]==1 && alpha->lconf[i2]==2 )
      {
        fr = GET_RAND;
        if(fr<.5) { //case 1,-1 -> -1,1
    
          p = swap(i1, i2);

          if( accept(p) ) {
            alpha->lconf[i1] = 2; alpha->lconf[i2] = 1;
            backup_data();
            wf *= p;
            accepted++; return true;
          }

        } else { //case 1,-1 -> 0,0

          p = crop(i1,i2);

          if( accept(p) ) {
            alpha->lconf[i1] = 0; alpha->lconf[i2] = 0; alpha->Na[0]+=2; alpha->Na[1]--; alpha->Na[2]--;
            backup_data();
            wf *= p;
            accepted++; return true;
          }

        }
        break;
      }

      if( alpha->lconf[i1]==2 && alpha->lconf[i2]==1 ) {
        fr = GET_RAND;
        if(fr<.5) { //case -1,1 -> 1,-1

          p = swap(i1, i2);

          if( accept(p) ) {
            alpha->lconf[i1] = 1; alpha->lconf[i2] = 2;
            backup_data();
            wf *= p;
            accepted++; return true;
          }

          alpha->lconf[i1] = 1; alpha->lconf[i2] = 2;
        } else { //case -1,1 -> 0,0

          p = crop(i1,i2);

          if( accept(p) ) {
            alpha->lconf[i1] = 0; alpha->lconf[i2] = 0; alpha->Na[0]+=2; alpha->Na[1]--; alpha->Na[2]--;
            backup_data();
            wf *= p;
            accepted++; return true;
          }

        }
        break;
      }
      // cases 1,0; -1,0; 0,1; 0,-1
      //just exchange the states

      //p = (d[alpha->lconf[i2]][i1]*d[alpha->lconf[i1]][i2])/(d[alpha->lconf[i1]][i1]*d[alpha->lconf[i2]][i2]);
      p = swap(i1, i2);

      if( accept(p) ) {
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        backup_data();
        wf *= p;
        accepted++; return true;
      }

      break;
    }
  } while (true);

  return false;
}


