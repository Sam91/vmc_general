#include "twoflavor.h"
#include "matrix.h"

/*
 * This virtual class is used to define a fully paired or unpaired state of two flavors (N0=N1) on any lattice.
 *
 * You need to derive a class that defines the matrix adx
 *
 */

twoflavor::twoflavor(int n, int l) : wavefunction( n, l )
{
  if(NS != 2) {
    cout << "ERROR: we can only do two flavors with this wave function\n";
    exit(-1);
  }  

  if(N%NS != 0) {
    cout << "N (number of sites) must be a multiple of " << NS << "\n";
    exit(-1);
  }

  N2 = N/2;

#if WFC
  adx = new complex<double>[N*N];
  current_inv = new complex<double>[N2*N2];
  old_inv = new complex<double>[N2*N2];
  sample_row = new complex<double>[N2];
#else
  adx = new double[N*N];
  current_inv = new double[N2*N2];
  old_inv = new double[N2*N2];
  sample_row = new double[N2];
#endif

  current_x = createint(NS, N2 );
  old_x = createint(NS, N2 ); 

  cff = 1.; //set the default wf normalization to 1
}

twoflavor::~twoflavor()
{
  cout << "twoflavor::~twoflavor()\n";
  delete[] adx;
  delete[] current_inv;
  delete[] old_inv;
  delete[] sample_row;

  destroy(current_x, NS);
  destroy(old_x, NS);
}

void twoflavor::backup_data()
{
  memcpy( old_x[0], current_x[0], N2*sizeof(int) );
  memcpy( old_x[1], current_x[1], N2*sizeof(int) );

#if WFC
  memcpy( old_inv, current_inv, N2*N2*sizeof(complex<double>) );
#else
  memcpy( old_inv, current_inv, N2*N2*sizeof(double) );
#endif

  wf_old = wf;
}

void twoflavor::restore_data()
{
  int* tmp;
#if WFC
  complex<double> *temp;
#else
  double* temp;
#endif

  tmp = current_x[0]; current_x[0] = old_x[0]; old_x[0] = tmp;
  tmp = current_x[1]; current_x[1] = old_x[1]; old_x[1] = tmp;

  temp = current_inv; current_inv = old_inv; old_inv = temp;

  wf = wf_old;
}

//get the wf = <alpha|psi>
void twoflavor::getwf()
{
  wf = 1.;

  int k0, k1;
  bool sgn = true;
  bool aux_sgn;

  //here, we have a single determinant and two spinon flavors, no holes

  k0 = k1 = 0; aux_sgn=true;
  for(int j=0; j<N; j++)
  {
    if( alpha->lconf[j]==0 )
    {
      current_x[0][k0] = j;
      k0++; aux_sgn = !aux_sgn;
    } else {
      current_x[1][k1] = j;
      k1++; if( aux_sgn==false ) sgn = !sgn;
    }
  }
  
  if( (k0!=N2) || (k1!=N2) ) {
    cout << "ERROR: incorrect lconf...\n";
    exit(-1);
  }

  //construct the d_matrix; we associate flavor 0 with the first (column), and flavor 1 with the second(row) index of the matrix
  for(int i=0; i<N2; i++)
    for(int j=0; j<N2; j++)
      current_inv[i*N2+j] = adx[ N*current_x[0][i] + current_x[1][j] ];

  //cout << "d_matrix :\n";
  //write_m(d_matrix, N2);

  wf = inverse(current_inv, N2 );

  if( abs(wf)>1e300 ) {
    cout << "ERROR(tf): divergent wf: " << wf << "\n";
    throw( -1 );
  }
  //cout << "sgn: " << sgn << "\n";
  if( !sgn ) wf = -wf;

  //cout << "wf norm: " << wf << "\n";

#if JASTROW
  wf *= jastrow();
#endif

  cout << "twoflavor::getwf(): " << std::scientific << wf << "\n";

  //return wf;
}

//swap operator of two sites i1 and i2 (corresponding to virtual_replacement)
//returns the ratio wf_after/wf_before
#if WFC
complex<double> twoflavor::swap(int i1, int i2, bool ratio)
{
  complex<double> r;
#else

double twoflavor::swap(int i1, int i2, bool ratio)
{
  double r;
#endif

  //cout << "Swapping " << i1 << ", " << i2 << " in\n";
  //alpha->print();

  //cout << "current_x:\n";
  //for(int i=0; i<N2; i++) cout << current_x[0][i] << ", " << current_x[1][i] << "\n";

  if( alpha->lconf[i1] == alpha->lconf[i2] )
  {
    //cout << "identical states in swap\n";
    return 1.;
  } 
    else
  {
    //if the current wf is a node, we cannot calculate the ratio.
    //In this case, we calculate from scratch and return 1.
    if( abs(wf)<SMALL )
    {
      if( ratio ) { cout << "ERROR: incorrect ratio for vanishing wf!\n"; exit(-1);}

      cout << "vanishing wf... restarting from scratch\n";
      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      r = wf;
      getwf(); //uses the current alpha configuration (needs to be swapped at that point)
      return wf/r;
    }
/*
    if( ratio ) {
      r = wf;
      backup_data();
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      getwf();
      c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      r = wf/r;
      restore_data();
      return r;
    } else
    {
      r = wf;
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      getwf();
      return wf/r;
    }
*/
    //re-arange such that i1=0, i2=1
    if( alpha->lconf[i1] == 1 ) {
      int ri = i1; i1 = i2; i2 = ri;
    }

    int row;

    //find the column number to update in the det
    for(row=0; (row<N2)&&(current_x[0][row]!=i1); row++);

    if( row>=N2 ) {
      cout << "Internal error i1 in swap\n";
      exit(-1);
    }

//    double **m = createdouble(N2);
//    double **m1 = createdouble(N2);

    for(int i=0; i<N2; i++) sample_row[i] = adx[ i2*N + current_x[1][i] ];

    if( ratio ) //here, we need to backup-restore since it is a single det...
    {
      backup_data();
      r = update_row( sample_row, row, current_inv, N2 );

      if( abs(r)< SMALL ) //intermediate det vanishes
      {
        cout << "Intermediate det vanishing (ratio)...\n";
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        r = wf_old;
        getwf();
        r = wf/r;
        restore_data();
        c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        return r;
      }
      current_x[0][row] = i2;
    } 
    else
    {

//      double det1,det2;
//      det1=inverse_mkeep( current_inv, m1, N2);
      //cout << "(1) det_before: " << 1./det1 << "\n";

      r = update_row( sample_row, row, current_inv, N2 );
      current_x[0][row] = i2;

/*      if(abs(r)>SMALL) {
        det2=inverse_mkeep( current_inv, m, N2 );
        //cout << "; det_after: " << 1./det2;
        //cout << "; r= " << r << "\n";
 
        if( (abs(det2)>SMALL && abs((det2*r - det1)/det2)>1e-5 ) ) {
          cout << "warn 1: d_matrix - old one " << alpha->lconf[i2] << " after swap of row " << row << " :\n";

          for(int i=0; i<N2; i++)
            for(int j=0; j<N2; j++)
              m[i][j] -= m1[i][j];
          write_m(m, N2);

          cout << "Replacement row:\n";
          for(int i=0; i<N2; i++) cout << adx[i2][current_x[1][i]]-adx[i1][current_x[1][i]] << ", ";
          cout << "\n";
        }
      }
*/
      if( abs(r)< SMALL ) { //intermetiate det vanishes
        cout << "Intermediate det vanishing...\n";
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        r = wf; 
        getwf();
        return wf/r;
      }
    }

    //find the row number to update in the det
    for(row=0; (row<N2)&&(current_x[1][row]!=i2); row++);
    
    if( row>=N2 ) {
      cout << "Internal error i2 in swap\n";
      exit(-1);
    }

    for(int i=0; i<N2; i++) sample_row[i] = adx[ N*current_x[0][i] + i1 ];

    if( ratio )
    {
      r *= det_ratio_column( sample_row, row, current_inv, N2 );
      
      restore_data();
    }
    else
    {
//      double det1,det2,r2;
//      det1=inverse_mkeep( current_inv, m1, N2 );
      //cout << "(2) det_before: " << 1./det1;

      r *= update_column( sample_row, row, current_inv, N2 );
//      r *= r2;
      current_x[1][row] = i1;

/*      if(abs(r2)>SMALL) {
        det2=inverse_mkeep( current_inv, m, N2 );
        //cout << "; det_after: " << 1./det2;
        //cout << "; r= " << r2 << "\n";
        //inverse_mkeep(current_inv, m, N2);

        if( (abs(det1)>SMALL) && abs((det2*r2 - det1)/det1)>1e-5 ) {
          cout << "warn 2: d_matrix-old one " << alpha->lconf[i1] << " after swap of row " << row << " :\n";
          for(int i=0; i<N2; i++)
            for(int j=0; j<N2; j++)
              m[i][j] -= m1[i][j];
          write_m(m, N2);

          cout << "Replacement row:\n"; 
          for(int i=0; i<N2; i++) cout << adx[current_x[0][i]][i1]-adx[current_x[0][i]][i2] << ", ";
          cout << "\n";
        }
      }
*/
      wf *= -r; //update the wf

      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
    }
//    destroy(m, N2 );
//    destroy(m1, N2 );
  }

  //if( abs(r-1.) > .01 )
  //  cout << "swap (" << alpha->lconf[i1] << ", " << alpha->lconf[i2] << ") returning " << r << "\n";

/*  if( !ratio ) {
    double wfbkp=wf;
    getwf();
    if( abs(1.-wfbkp/wf)>1e-4 ) {
      cout << "ERROR: wrong wf update: wf=" << wf << "; upd wf=" << wfbkp << "\n";
    }
  } //if conf1 == conf2 loop
*/
  return -r;
}

void twoflavor::find_starting_conf()
{
  int* Na0 = new int[NS];

  for(int n=0; n<NS; n++) Na0[n] = N/NS;

  alpha->set_random_conf( Na0 );

  while( true )
  {
    try {
      getwf();
    } catch( int e )
    {
      cout << "Decreasing normalization...\n";
      correct_cff( false );
      continue;
    }
    break;
  }

  //cout << "Max config:\n";
  //alpha->print();
  cout << "wf: " << wf << "\n";

  delete[] Na0;
}

#if WFC
complex<double> twoflavor::crop(int i1,int i2,bool b) {return 0.;}
#else
double twoflavor::crop(int i1,int i2,bool b) {return 0.;}
#endif

//corrects the wave fuction by pow(10,s 50), depending on the sign of s
void twoflavor::correct_cff(bool s)
{
  double r;

  if( s )
    r = pow(10., +40./(double)N2 );
  else
    r = pow(10., -60./(double)N2 );

  normalize( r );

  cout << "TF: cff changed to " << cff << endl;
  //create_ad();
}

void twoflavor::normalize(double r)
{
  cout << "Normalizing with r = " << r << endl;

  for(int j=0; j<N*N; j++) adx[j] *= r;
  cff *= r;
}

