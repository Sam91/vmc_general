/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of first, second, and third neighbor correlator without averaging over rotations (works for kagome lattice)
 * Note that we average over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = 15
 */


  for(int nx=0; nx<L; nx++ ) //loop over carthesians
  {
    for(int ny=0; ny<L; ny++ )
    {
      for(int q=0; q<Q; q++)
      {
        int i1 = alpha->mylattice->j(nx,ny,q);

        //swap with ratio=true(i.e., only the ratio is calculated but the wf is not updated)

        //nearest neighbors
#if WFC
        fj[2*q  ] += swap( i1, alpha->mylattice->links[i1][0][1], true ).real();
        fj[2*q+1] += swap( i1, alpha->mylattice->links[i1][0][2], true ).real();
#else
        fj[2*q  ] += swap( i1, alpha->mylattice->links[i1][0][1], true );
        fj[2*q+1] += swap( i1, alpha->mylattice->links[i1][0][2], true );
#endif

        //next neighbors
#if WFC
        fj[2*Q+2*q  ] += swap( i1, alpha->mylattice->links[i1][1][1], true ).real();
        fj[2*Q+2*q+1] += swap( i1, alpha->mylattice->links[i1][1][2], true ).real();
#else
        fj[2*Q+2*q  ] += swap( i1, alpha->mylattice->links[i1][1][1], true );
        fj[2*Q+2*q+1] += swap( i1, alpha->mylattice->links[i1][1][2], true );
#endif

        //third neighbors
#if WFC
        fj[4*Q+q] += swap( i1, alpha->mylattice->links[i1][2][1], true ).real();
#else
        fj[4*Q+q] += swap( i1, alpha->mylattice->links[i1][2][1], true );
#endif

      }
    }
  } // end loop over all sites

