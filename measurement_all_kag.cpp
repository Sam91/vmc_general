/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of the full correlator between all sites (first used for the Kagome lattice Dirac state)
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = Q*Q*L2
 */

  //this should be fixed to use linear indices directly

  int xx2, yy2, j1, j2;

  for (int x1=0; x1<L; x1++)
  {
    for (int y1=0; y1<L; y1++)
    {
      for (int q1=0; q1<Q; q1++)
      {
        j1 = q1 + Q*( L*x1 + y1);

        for (int x2=0; x2<L; x2++)
        {
          xx2 = x1+x2; 
          if( xx2>L-1 ) xx2 -= L;

          for (int y2=0; y2<L; y2++)
          {
            yy2 = y1+y2;
            if( yy2>L-1 ) yy2 -= L;

            for (int q2=0; q2<Q; q2++)
            {
              j2 = q2 + Q*( L*xx2 + yy2);

              if( j1==j2 )
                fj[ q1 ] += .75;
              else
                fj[ q1 + Q*(q2 + Q*(L*x2+y2)) ] += swap( j1, j2, true);

            }
          }
        }
      }
    }
  }

