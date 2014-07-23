/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of sum over sublattice corelators
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = 3
 */


  for(int i1=0; i1<LD; i1++ ) //loop over all pairs of sublattices
  {
    for(int i2=0; i2<LD; i2++ )
    {
      for(int q=0; q<Q; q++)
      {
        if( i1==i2 ) {
          fj[q] += 2.;
          continue;
        }
#if WFC
        fj[q] += swap( Q*i1 + q, Q*i2 + q, true ).real();  //swap with ratio=true(i.e., only the ratio is calculated but the wf is not updated)
#else
        fj[q] += swap( Q*i1 + q, Q*i2 + q, true );
#endif
      }
    }  
  } // end loop over all sites

