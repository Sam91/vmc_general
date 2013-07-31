/*
 * Measurement that can be included in wavefunction::accumulate()
 *
 * Implements measurement of first, second, and third neighbor correlator (first used for triangular lattice states)
 * Also, 4-site ring terms are measured (defined in lattice.cpp)
 * Note that we are averaging over lattice translations here, so wavefunction::average() needs to divide by L2
 *
 * For this measurement, we need NO = 5
 */


  for(int i1=0; i1<N; i1++ ) //loop over all sites
  {

    //nearest neighbors
/*
     for(int n=1; n<=alpha->mylattice->links[i1][0][0]; n++)
    {
#if WFC
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true ).real();  //swap with ratio=true(i.e., only the ratio is calculated but the wf is not updated)
#else
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][n], true );
#endif
    }
*/
    //
    //anisotropic nearest neighbor (triangular lattice only)
#if WFC
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][1], true ).real();
      fj[4] += swap( i1, alpha->mylattice->links[i1][0][2], true ).real();
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][3], true ).real();
#else
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][1], true );
      fj[4] += swap( i1, alpha->mylattice->links[i1][0][2], true );
      fj[0] += swap( i1, alpha->mylattice->links[i1][0][3], true );
#endif

    //next-nearest neighbor exchange
    for(int n=1; n<=alpha->mylattice->links[i1][1][0]; n++)
    {
#if WFC
      fj[1] += swap( i1, alpha->mylattice->links[i1][1][n], true ).real();
#else
      fj[1] += swap( i1, alpha->mylattice->links[i1][1][n], true );
#endif
    }

    //third-neighbor exchange
    for(int n=1; n<=alpha->mylattice->links[i1][2][0]; n++)
    {
#if WFC
      fj[2] += swap( i1, alpha->mylattice->links[i1][2][n], true ).real();
#else
      fj[2] += swap( i1, alpha->mylattice->links[i1][2][n], true );
#endif
    }

#if (NS==2) && (NO>2)

    //four-site plaquette term (this measurement is optimized for TWO flavors, where a four-cycle does not arise)
    //loop over plaquettes
//cout << "Looping " << alpha->mylattice->plaquette4[i1][0][0] << " plaquettes.\n";
    for(int n=1; n<=alpha->mylattice->plaquette4[i1][0][0]; n++)
    {
//cout << "n: " << n << "\n";

      i2 = alpha->mylattice->plaquette4[i1][0][3*n-2];
      i3 = alpha->mylattice->plaquette4[i1][0][3*n-1];
      i4 = alpha->mylattice->plaquette4[i1][0][3*n-0];

      //cout << "Plaquette: " << "(" << i1 << "," << i2 << "," << i3 << "," << i4 << ")\n";

      if( (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(aaaa)
      {
        fj[3] += 1.;
        continue;
      }

      if( (alpha->lconf[i1]!=alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(abbb)
      {
//        r = swap(i1, i2, true);
//        if(abs(r)>50.) cout << "ERROR r1a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i1, i2, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i1, i2, true);
#endif
        continue;
      }

      if( (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) ) //(babb)
      {
//        r = swap(i2, i3, true);
//        if(abs(r)>50.) cout << "ERROR r1b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i2, i3, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i2, i3, true);
#endif
        continue;
      }

      if( (alpha->lconf[i3]!=alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) && (alpha->lconf[i1]==alpha->lconf[i2]) ) //(bbab)
      {
//        r = swap(i3, i4, true);
//        if(abs(r)>50.) cout << "ERROR r1c, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i3, i4, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i3, i4, true);
#endif
        continue;
      }

      if( (alpha->lconf[i4]!=alpha->lconf[i1]) && (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]==alpha->lconf[i3]) ) //(bbba)
      {
//        r = swap(i4, i1, true);
//        if(abs(r)>50.) cout << "ERROR r1d, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i4, i1, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i4, i1, true);
#endif
        continue;
      }

      if( (alpha->lconf[i1]==alpha->lconf[i2]) && (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]==alpha->lconf[i4]) ) //(aabb)
      {
//        r = swap(i1, i3, true);
//        if(abs(r)>50.) cout << "ERROR r2a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i1, i3, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i1, i3, true);
#endif
        continue;
      }

      if( (alpha->lconf[i2]==alpha->lconf[i3]) && (alpha->lconf[i3]!=alpha->lconf[i4]) && (alpha->lconf[i4]==alpha->lconf[i1]) ) //(abba)
      {
//        r = swap(i2, i4, true);
//        if(abs(r)>50.) cout << "ERROR r2b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//        fj[2] += r.real()
        fj[3] += swap(i2, i4, true).real();
#else
//        fj[2] += r;
        fj[3] += swap(i2, i4, true);
#endif
        continue;
      }

      if( (alpha->lconf[i1]!=alpha->lconf[i2]) && (alpha->lconf[i2]!=alpha->lconf[i3]) && (alpha->lconf[i3]!=alpha->lconf[i4]) ) //(abab); this is a product of two transpositions
      {
        backup_data();
        int c;

        r = swap(i1, i2, false); //here, we update the matrix
//        if(abs(r)>50.) cout << "ERROR r4a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";

        if( abs(r)<SMALL ) //intermediate wf vanishes
        {
          c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
          r = wf_old;
          getwf();
//          r = wf/r;
//          if(abs(r)>50.) cout << "ERROR r4a, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";
#if WFC
//          fj[2] += r.real();
          fj[3] += (wf/r).real();
#else
//          fj[2] += r;
          fj[3] += wf/r;
#endif
          restore_data();
          c = alpha->lconf[i1]; alpha->lconf[i1]=alpha->lconf[i2]; alpha->lconf[i2]=c;
          c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
          continue;
        }

        r *= swap(i3, i4, false); //here, we update the matrix
//        if(abs(r)>50.) cout << "ERROR r4b, " << n << ", (" << i1 << "," << i2 << "," << i3 << "," << i4 << "): " << r << "\n";

#if WFC
        fj[3] += r.real();
#else
        fj[3] += r;
#endif
        restore_data();
        c = alpha->lconf[i1]; alpha->lconf[i1]=alpha->lconf[i2]; alpha->lconf[i2]=c;
        c = alpha->lconf[i3]; alpha->lconf[i3]=alpha->lconf[i4]; alpha->lconf[i4]=c;
      }
    } //plaquette loop
#endif //NS==2

  } // end loop over all sites

