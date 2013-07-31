/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of first neighbor correlator
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = 1
 */


  for(int i1=0; i1<N; i1++ ) //loop over all sites
  {

    //nearest neighbors
     for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)
    {
#if WFC
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true ).real();  //swap with ratio=true(i.e., only the ratio is calculated but the wf is not updated)
#else
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true );
#endif
    }

  } // end loop over all sites

