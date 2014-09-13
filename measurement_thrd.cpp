/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of all three third neighbor correlators (kagome)
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = 3
 */


  for(int i1=0; i1<N; i1++ ) //loop over all sites
  {
//     int q = i1%3;

    //first neighbor
    for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)
    {
#if WFC
      fj[ n-1 + 2*i1 ] += swap( i1, alpha->mylattice->links[i1][0][n], true ).real();
#else
      fj[ n-1 + 2*i1 ] += swap( i1, alpha->mylattice->links[i1][0][n], true );
#endif
    }


    //third-neighbor
    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)
    {
#if WFC
      fj[ 2*N + i1 ] += swap( i1, alpha->mylattice->links[i1][2][n], true ).real();
#else
      fj[ 2*N + i1 ] += swap( i1, alpha->mylattice->links[i1][2][n], true );
#endif
    }

  } // end loop over all sites

