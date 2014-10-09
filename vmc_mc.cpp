/***************************************************************************
 *   Copyright (C) 2011 by Samuel Bieri                                    *
 *   samuel.bieri@a3.epfl.ch                                               *
 *                                                                         *
 ***************************************************************************/
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "vmc.h"

/* !!!!!!!!!
 * THIS CODE IS DEPRICATED AND UNUSED! I keep it here just as a reminder...
 * !!!!!!!!!
 */

/* In this file, we define code related to Monte Carlo walking in the space of Ising states */ 

void vmc::randomize()
{
  //srand( time(0) );
  long sek;
  time(&sek);
  srand( (unsigned)sek );
  //srand( 0 );
}

void vmc::set_random_conf()
{
  int h1, s1, s2;
  double r,x1,x2;
  int ii,jj;

  //cout << "starting new set_random_conf...\n";

  /* Set a random configuration */
  h1=LD-SPINS_UP-SPINS_DN;
  s1=SPINS_UP;
  s2=SPINS_DN;
  //cout << " set random conf\n";

  for(ii=0;ii<L;ii++) {
    for(jj=0;jj<L;jj++) {
      r  = GET_RAND;
      x1 = (double)h1 / (double)(h1+s1+s2);
      x2 = (double)(h1+s1) / (double)(h1+s1+s2);
      //cout << "x1 x2: " << x1 << " " << x2 << "\n";
      //printf( "%i, %i: ", ii, jj);
      if (r<x1) {
        conf[ii][jj] = 0; h1--;
      } else {
        if (r<x2) {
          conf[ii][jj] = 1;
          s1--;
        } else {
          conf[ii][jj] =-1;
          s2--;
        }
      }
      //cout << conf[ii][jj] << " ";
    }
    //cout << "\n";
  }
  ZERO_STATUS = false;

/*  for(ii=0; ii<L; ii++)
    for(jj=0; jj<L; jj++)
      conf[ii][jj] = ((ii+2*jj)%3)-1;
*/
  getwf();
  //cout << "set_random_conf: wf=" << wf << " ("<< det1 << ", "<< det3 <<")\n";
}

//Set a new random spin configuration in conf
int vmc::find_random_conf()
{
  //cout << "starting new find_random_conf...\n";

  /* Try to get a non-vanishing configuration */
  for(int n=0; n<LD; n++)
  {
    set_random_conf();
    thermalize();
#if XYODD
    cout << "find_random_conf: wf=" << std::scientific << setprecision(2) << wf << " ("<< det1 << ", "<< det2 << ", " << det3 <<")\n";
#else
#if ONESLATTERDET
    cout << "find_random_conf: wf=" << std::scientific << setprecision(2) << wf << "\n";
#endif
#if TWOSLATTERDETS || XYDET
    cout << "find_random_conf: wf=" << std::scientific << setprecision(2) << wf << " ("<< det1 << ", "<< det3 <<")\n";
#endif
#endif

    if( !ZERO_STATUS && pow(abs(wf),2)>ALMOST_ZERO) {
      //cout << "find_random_conf: found suitable starting config:\n";
      //write_m(conf, L);
      return 0;
    }
  }

  //cout << "d3_matrix:\n";
  //write_m(d3_matrix,H);

  //cout << "inv_m3:\n";
  //write_m(current_inv3,H);

  //No conf found with non-zero overlapp (the wf possibly vanishes)
  return -1;
}

//Set a given ising configuration -- used for testing purpose (or exact projection)
void vmc::set_conf( int **c )
{
  //cout << "setting conf\n";
  for(int i=0;i<L;i++) {
    for(int j=0;j<L;j++) {
      //cout << "("<< i << ","<< j<< ")";
      conf[i][j] = c[i][j];
    }
    //cout << "\n";
  }
  getwf();
}

/* experience shows that the system is equilibrated (wf^2 saturates) after roughly twice systemSize steps*/
void vmc::thermalize()
{
  int nn = 4*N;
  cout << "thermalizing with " << nn << " steps\n";
  for(int i=0; i<nn; i++ )
    walk_one_step();
}

/* We speak of a thermalized monte-carlo dynamics, as soon as the wf we sample has saturated to some large norm. 
 * To visualize this, we can plot ln(|wf|^2) as a function of MC steps to see approximately when it saturates.
 * Another important quantitiy is the acceptance rate. Of course we would like it to be as high as possible
 * for fast MC ensemble generation. Unfortunately, it seems to go down to 25-30%.
*/
void vmc::visualize_thermalization()
{
  ofstream f_ar("ar.txt");
  ofstream f_wf("wf.txt");

  ACCEPTED_MOVES=0;
  STEPS_TOTAL=0;
  for(int j=0; j<1000; j++) { // perform a MC run and get the data
    //cout << " " << log( wf * wf );
    walk_one_step();
    f_ar << ar() << "\n";
    f_wf << log( wf * wf ) << "\n";
  }

  f_ar.close();
  f_wf.close();
}

/* a walk has number of steps of order of system size (2*L*L)
 * Approximately every electron has been exchanged once.
 */
void vmc::run_a_walk()
{
  //cout << "running a vmc walk\n";
  STEPS_TOTAL = 0;
  ACCEPTED_MOVES = 0;
  for( int i=0;i<4*LD;i++ ) walk_one_step();

  //cout << "MC steps: " << STEPS_TOTAL << "\n";
//  cout << "Accepted moves: " << ACCEPTED_MOVES << "/" << STEPS_TOTAL << "\n";
//  write_m(conf,L);

  if(ACCEPTED_MOVES==0)
  {
    cout << "WARN: we are stuck in the following config:\n";
    write_m(conf,L);
  }

  if( ACCEPTED_MOVES==STEPS_TOTAL )
  {
    cout << "WARN: all steps accepted... There is probably an issue somewhere.\n";
    cout << "wf = " << wf << "; old_wf = " << old_wf[0] << "\n";
    cout << "det1 = " << det1 << "; old_det1 = " << old_det1[0] << "\n";
    cout << "det3 = " << det3 << "; old_det3 = " << old_det3[0] << "\n";
    exit(-1);
  }

}

/* fraction of accepted MC moves */
double vmc::ar() { return (double)ACCEPTED_MOVES / (double)STEPS_TOTAL; }

//this function is deprecated; we do this now in wavefunction.cpp
void vmc::walk_one_step()
{
  STEPS_TOTAL++;

  // choose at random two neighboring points on the lattice to exchange.
  //fr = (double)::std ::rand()/RAND_MAX;
  //fr = rg.Random();
  fr = GET_RAND;

  ir = (int) (fr*LD);

  i1 = (int) (ir/L);                /* i1,j1 is the first point to exchange */
  j1 =       (ir%L);
  //Unfortunately, this probabilty is inclusive 1, so we have to fix this case. Like that, we slightly overweight border-hopping.
  //if(i1>=L) {i1=L-1; cout << "inclusiv problem\n";}
  //if( fr > 1.0 || fr < 0.0 ) cout << "fr: " << fr << "\n";
  
  //fr = (double)::std::rand()/RAND_MAX;
  //fr = rg.Random();

// pick another random point until we find a different flavor
  do {
    fr = GET_RAND;
  
//Choose a neares neighbor at random (local updates only)
/*  if( fr<0.33333333 ) { // i2,j2 is the second point to exchange
    i2=i1; j2=j1+1; if(j2==L) j2=0;
  } else
    if( fr<0.66666666 ) {
      j2=j1; i2=i1+1; if(i2==L) i2=0;
    } else {
      i2=i1+1; j2=j1+1; if(i2==L) i2=0; if(j2==L) j2=0;
    }
*/
    //Randomy chose a second point (nonlocal updates)
    ir = (int) (fr*LD);
    i2 = (int) (ir/L);			// i2,j2 is second point to exchange
    j2 =       (ir%L);

    //Print out the two points to swap
    //cout << "x1=("<< i1 <<","<< j1 <<"), x2=("<< i2 << ","<< j2<<").\n";
 
  } while( conf[i1][j1] == conf[i2][j2] );

  /* Continue if identical flavors occupy the sites (in this case the wf remains unchanged) */
  //if( conf[i1][j1] == conf[i2][j2] ) return;

  //cout << "x1=("<< i1 <<","<< j1 <<"), x2=("<< i2 << ","<< j2<<"): "<< conf[i1][j1] << ", "<< conf[i2][j2] << ".\n";
  /*
  cout << "conf:\n";
  write_m(conf,L);
  

  // print out the nz position
  for(int i=0; i<H; i++)
    cout << "("<< current_nz_x[i] << ","<<current_nz_y[i] << ")";
  cout << "\n";
  */

  // backup the data because we are going to perform a move
  backup_data(); // sets old_wf = wf
  virtual_replacement(i1,j1,i2,j2); //Updates the wf

#if DEBUG
  //check if we get the same wf as from scratch
  complex<double> wf_test = wf; 
#if DET1C
  complex<double> det1_test = det1; 
#else
  double det1_test = det1; 
#endif
#if DET3C
  complex<double> det3_test = det3; 
#else
  double det3_test = det3; 
#endif
  getwf();

  if( abs((wf-wf_test)/wf) > 1e-2 && (abs(wf) > ALMOST_ZERO) && (abs(wf_test) > ALMOST_ZERO) && !ZERO_STATUS) {
    cout << "gss, we have a problem: " << wf << " " << wf_test << " " << abs((wf-wf_test)/wf) << "; x1=(" << i1 << "," << j1 << "), x2=(" << i2 << ","<< j2 << "): "<< conf[i1][j1] << ", "<< conf[i2][j2] << ".\n";
    cout << "det1: " << det1 << " " << det1_test << "\n";
    cout << "det3: " << det3 << " " << det3_test << "\n";
    cout << "ZERO_STATUS: " << ZERO_STATUS << "\n";

    //cout << "old_inv3:\n";
    //write_m(old_inv3[0], H);

    //cout << "current_inv3:\n";
    //write_m(current_inv3, H);

  } /*else { 
    cout << "Here it matches...\n";
    cout << wf_test << ", " << wf << "\n";
  }*/
#endif

  /* Check also the two determiniants; note, however, that we do not get the sign right here. But we only need it for measurment purpose.
  if( (!ZERO_STATUS) && abs(det1_test-det1) > 1e-10 ) cout << "det1 mismatch: " << det1_test << "; " << det1 << "\n";
  if( (!ZERO_STATUS) && abs(det3_test-det3) > 1e-10 ) cout << "det3 mismatch: " << det3_test << "; " << det3 << "\n";
  */

  //check if inv_adx3 is still skew symmetric
/*  if( !checkskew(current_inv3, H) ) {
    cout << 
  }
*/
  /* check if the wf is given by the product of the dets */
/*  if(ZERO_STATUS) cout << "ZERO_STATUS is true\n";
  if(sign) { 
    if( abs(wf - det1*det3)>1e-10 ) cout << "Internal error with inconsistent dets (+)... check this! wf=" << wf << "; det1*det3="<< det1*det3 << "\n";
  } else
    if( abs(wf + det1*det3)>1e-10 ) cout << "Internal error with inconsistent dets (-)... check this! wf=" << wf << "; det1*det3="<< det1*det3 << "\n";
*/
  if( INF_STATUS ) //if we get a very large value for the new wf, then the current wf is probably already wrong. Restart from scratch
  {
    v=conf[i2][j2]; conf[i2][j2]=conf[i1][j1]; conf[i1][j1]=v;
    restore_data();

    getwf(); //recalculate and start from scratch
    INF_STATUS = false;
    backup_data();
    virtual_replacement(i1,j1,i2,j2);
  }

  if( ZERO_STATUS ) // Reject nodes immediately
  {
    v=conf[i2][j2]; conf[i2][j2]=conf[i1][j1]; conf[i1][j1]=v;
    restore_data();
    ZERO_STATUS = false;
    return;
  }

#if JASTROW
  p  = abs(wf) / abs(old_wf[0]) * jastrow(i1,j1,i2,j2);
#else
  p  = abs(wf) / abs(old_wf[0]);
#endif

  //cout << "p: " << p << "\n";

  p *= p;

  if( p<1. )
  {
    //fr = (double)::std::rand()/RAND_MAX;
    //fr = rg.Random();
    fr = GET_RAND;
    if( fr>p ) {
      /* Do not perform exchange */
      v=conf[i2][j2];
      conf[i2][j2]=conf[i1][j1];
      conf[i1][j1]=v;
      restore_data();
      //cout << "rejected\n";
      return;
    }
  }

  /* Otherwise -- perform exchange */
  //cout << "accepted. wf= "<< wf <<"\n";
  ACCEPTED_MOVES++;
  return;
}

