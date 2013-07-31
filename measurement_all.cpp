/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of the full correlator between all sites (first used for the Kagome lattice Dirac state)
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = Q*Q*L2
 */

  /* simple nearest-neighbor S.S
  for (int i1=0; i1<N;i1++){
    for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)
    {
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true );
    }
  }
  */

  int i1r, i2r;
  for (int l1=0; l1<Q; l1++){ // loop over first label within unit cell
    for (int l2=0; l2<Q; l2++){ // over second label
      for (int v1=0; v1<L; v1++){ // loop over vector x-coordinate
        for (int v2=0; v2<L; v2++){ // loop over vector y-coordinate
          for (int i1=0; i1<L; i1++){ // loop over site x-coordinate
            for (int i2=0; i2<L; i2++){ // loop over site y-coordinate
              i1r = i1+v1;
              i2r = i2+v2;
              alpha->mylattice->torus(i1r, i2r);
              fj[ (Q*l2+l1)*L2+L*v2+v1 ] += swap(Q*(L*i2+i1)+l1, Q*(L*i2r+i1r)+l2, true);
            }
          }
        }
      }
    }
  }

