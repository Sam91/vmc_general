//this now depends on the particular hilbert space we want to consider

#include "isingstate.h"

//This generates the space na=3 with the constraint N1 = N2 and 00->12,21; 12->00,21; 21->00,12
void isingstate::step()
{
  int i1, i2;
  do {

    //pick two sites at random
    double fr = N*GET_RAND;

    i1 = (int)(fr);
    i2 = (int)((fr-i1)*N);

    if( i1==i2 ) continue;

    if( lconf[i1]==lconf[i2] )
    {
      if( lconf[i1]==0 ) {
/*        fr = GET_RAND;
        if(fr<.5) { lconf[i1] = 1; lconf[i2] = 2; Na[0]-=2; Na[1]++; Na[2]++;}
        else { lconf[i1] = 2; lconf[i2] = 1; Na[0]-=2; Na[1]++; Na[2]++; }
*/
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

