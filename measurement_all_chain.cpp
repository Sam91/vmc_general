/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of the full correlator between all sites (first used for the Kagome lattice Dirac state)
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = Q*Q*LD
 */

  for( int j=0; j<N; j++ )
  {
    for( int k=0; k<N-1; k++)
    {
      int j2 = j+k+1; if(j2>=N) j2 -= N;
      fj[ k ] += swap(j, j2, true);
    }
  }

