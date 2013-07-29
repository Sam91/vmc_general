#include "u1.h"
#include "matrix.h"

//define a determinental U(1) state on any lattice for NS flavors (wf is a product of NS determinants).

u1::u1(int n, int l) : wavefunction( n, l )
{
  
  if(N%NS != 0) {
    cout << "N (number of sites) must be a multiple of " << NS << "\n";
    exit(-1);
  }

  pars = new parameters;
  bpars = pars;
  pars->t = new (double(**[NS]));
  pars->NF = new int[NS];

  adx = new (double(*[NS]));
  current_inv = new (double(*[NS]));
  old_inv = new (double(*[NS]));

  current_x = createint(NS, N);
  old_x = createint(NS, N); 

  for(int n=0; n<NS; n++)
  {
    pars->t[n] = createdouble( N );
    pars->NF[n] = N/NS;
    adx[n] = new double[N*N];
    current_inv[n] = new double[N*N];
    old_inv[n] = new double[N*N];
  }
#if WFC
  wfn = new complex<double>[NS];
  wfn_old = new complex<double>[NS];
#else
  wfn = new double[NS];
  wfn_old = new double[NS];
#endif
  cff = new double[NS];
  for(int i=0; i<NS; i++) cff[i] = 1.;

  pars->t1 = 1.;
  pars->t2 = 1.;
  pars->ap = true;
  pars->N = N;
}

u1::~u1()
{
  for(int n=0; n<NS; n++)
  {
    destroy(pars->t[n], N);
    delete[] adx[n];
    delete[] current_inv[n];
    delete[] old_inv[n];
  }
  delete[] adx;
  delete[] current_inv;
  delete[] old_inv;

  delete[] pars->t;
  delete[] pars->NF;
  bpars = NULL;
  delete pars;

  delete[] cff;
  delete[] wfn;
  delete[] wfn_old;

  destroy(current_x, NS);
  destroy(old_x, NS);
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
void u1::set_hopping(double t, double tt)
{
  pars->t1 = t;
  pars->t2 = tt;

  int i2;

  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n=0; n<NS; n++)
        pars->t[n][i][j] = 0.;

  for(int i1=0; i1<N; i1++) {
    for(int n=0; n<NS; n++) {
      for(int k=1; k<=alpha->mylattice->connectivity[i1][0][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][0][k];
        if( pars->ap && ( abs(i1%L - i2%L)>L/2 ) ) pars->t[n][i1][i2] += -t; //FIXME: not sure if BC work for lattices with more than one site per unit cell
        else pars->t[n][i1][i2] += t;
      }
      for(int k=1; k<=alpha->mylattice->connectivity[i1][1][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][1][k];
        if( pars->ap && ( abs(i1%L - i2%L)>L/2 ) ) pars->t[n][i1][i2] += tt; //FIXME: same as above
        else pars->t[n][i1][i2] += -tt;
      }
    }
  }
}

void u1::print()
{
  cout << "wf pars->\n";
  cout << " J1, J2 = " << js[0] << ", " << js[1] << "\n";

  //cout << "t: " << t << "\n";
}

//initialize the U(1) state
void u1::create_ad()
{
  pars->desc = "U(1)";

  double **m = createdouble( N );
  double *v = new double[ N ];

  for(int n=0; n<NS; n++)
  {
    copy_m(m, pars->t[n], N);
    eigvects(m, v, N);
    cout << "v[" << n << "]: ";
    for(int i=0; i<N; i++) cout << v[i] << "; ";
    cout << "\n";

    //take the lowest Na eigenvectors
    int i;
    for(i=0; i<pars->NF[n]; i++)
      for(int j=0; j<N; j++) adx[n][j*pars->NF[n]+i] = m[i][j]; //note that the eigenvalue index is the last index in adx!

    if(i<N) {
      //cout << "abs diff " << abs(v[i]-v[i-1]) << "\n";
      if( abs(v[i]-v[i-1]) < 1e-8 )
        cout << "WARN: degenerate state\n";
      cout << "Fermi energy is " << v[i-1] << ", " << v[i] << " at position " << i-1 << "\n"; 
    }
  }

  destroy(m, N);
  delete[] v;
  normalize( 2. );
}

void u1::backup_data(int flavor1, int flavor2)
{
  memcpy( old_x[flavor1], current_x[flavor1], pars->NF[flavor1]*sizeof(int) );
  memcpy( old_x[flavor2], current_x[flavor2], pars->NF[flavor2]*sizeof(int) );

  memcpy( old_inv[flavor1], current_inv[flavor1], pars->NF[flavor1]*pars->NF[flavor1]*sizeof(double) );
  memcpy( old_inv[flavor2], current_inv[flavor2], pars->NF[flavor2]*pars->NF[flavor2]*sizeof(double) );

  bkp_f1 = flavor1;
  bkp_f2 = flavor2;

  wf_old = wf;
  wfn_old[flavor1] = wfn[flavor1];
  wfn_old[flavor2] = wfn[flavor2];
}

void u1::backup_data() //This backup is needed for the 4 site plaquette exchange term (which is only defined for TWO flavors)
{
  backup_data(0, 1);
}

void u1::restore_data()
{
  int* tmp;
  double* temp;

  tmp = current_x[bkp_f1]; current_x[bkp_f1] = old_x[bkp_f1]; old_x[bkp_f1] = tmp;
  tmp = current_x[bkp_f2]; current_x[bkp_f2] = old_x[bkp_f2]; old_x[bkp_f2] = tmp;

  temp = current_inv[bkp_f1]; current_inv[bkp_f1] = old_inv[bkp_f1]; old_inv[bkp_f1] = temp;
  temp = current_inv[bkp_f2]; current_inv[bkp_f2] = old_inv[bkp_f2]; old_inv[bkp_f2] = temp;

  wf = wf_old;
  wfn[bkp_f1] = wfn_old[bkp_f1];
  wfn[bkp_f2] = wfn_old[bkp_f2];
}

//get the wf = <alpha|psi>
void u1::getwf()
{
  wf = 1.;

  int k;
  bool sgn = true;
  bool aux_sgn;

  for(int n=0; n<NS; n++) //loop over flavors
  {
    k=0; aux_sgn=true;
    for(int j=0; j<N; j++) {//loop over all sites and pick the ones occupied by the flavor
      if( alpha->lconf[j] == n ) {
        for(int i=0; i<pars->NF[n]; i++) { //loop over eigenvectors
          current_inv[n][k*pars->NF[n]+i] = adx[n][j*pars->NF[n]+i];
        }
        current_x[n][k] = j;
        if( !aux_sgn ) {sgn = !sgn;}
        k++; //increment space index
      } else {
        if( alpha->lconf[j] > n ) aux_sgn = !aux_sgn; //permutation sign
      }
    }
    //cout << "d_matrix for " << n << ":\n";
    //write_m(current_inv[n], pars->NF[n]);

    wfn[n] = inverse(current_inv[n], pars->NF[n] );
    wf *= wfn[n];

    //cout << "inverse: \n";
    //write_m(current_inv[n], pars->NF[n]);
  }

  if( abs(wf)>1e300 ) {
    cout << "Warn: divergent wf\n";
  }
  //cout << "sgn: " << sgn << "\n";
  if( !sgn ) wf = -wf;

/*  for(int n=0; n<NS; n++) {
    cout << "Flavor " << n << ": ";
    for(int i=0; i<pars->NF[n]; i++) cout << current_x[n][i] << ", ";
    cout << "\n";
  }
*/
  //cout << "wf norm: " << wf << "\n";

#if JASTROW
  wf *= jastrow();
#endif

  cout << "getwf(): " << std::scientific << wf << "; ";
  for(int n=0; n<NS; n++) cout << wfn[n] << ", ";
  cout << endl;
  
  //return wf;
}

//swap operator of two sites i1 and i2 (corresponding to virtual_replacement)
#if WFC
complex<double> u1::swap(int i1, int i2, bool ratio)
{
  complex<double> r;
#else
double u1::swap(int i1, int i2, bool ratio)
{
  double r;
#endif

  if( alpha->lconf[i1] == alpha->lconf[i2] ) {
    //cout << "identical states in swap\n";
    return 1.;
  } 
    else
  {

    //if the current wf is a node, we cannot calculate the ratio. In this case, we calculate from scratch and return 1.
    if( abs(wf)<1e-200 )
    {
      if( ratio ) { cout << "ERROR: incorrect ratio for vanishing wf!\n"; exit(-1);}

      cout << "vanishing wf... restarting from scratch\n";
      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      getwf(); //uses the current alpha configuration (needs to be swapped at that point)
      return 1.;
    }
/*
    cout << "swap(" << i1 << ", " << i2 << ", " << ratio << ")\n";

    cout << "f-" << alpha->lconf[i1] << ": ";
    for(int i=0; i<pars->NF[alpha->lconf[i1]]; i++) cout << current_x[alpha->lconf[i1]][i] << ", ";
    cout << "\n";

    cout << "f-" << alpha->lconf[i2] << ": ";
    for(int i=0; i<pars->NF[alpha->lconf[i2]]; i++) cout << current_x[alpha->lconf[i2]][i] << ", ";
    cout << "\n";
*/
    int row;

    //find the row number of the first det
    for(row=0; (row<pars->NF[alpha->lconf[i1]])&&(current_x[alpha->lconf[i1]][row]!=i1); row++);

    if( row>=pars->NF[ alpha->lconf[i1] ] ) {
      cout << "Internal error i1 in swap\n";
      exit(-1);
    }

    //double **m = createdouble(pars->NF[0]);
    //double **m1 = createdouble(pars->NF[0]);

    if( ratio )
    {
      r = det_ratio_row( &(adx[alpha->lconf[i1]][i2*pars->NF[alpha->lconf[i1]]]), row, current_inv[alpha->lconf[i1]], pars->NF[alpha->lconf[i1]] );

      if( abs(r)< SMALL ) return 0.;
    } 
    else
    {

/*      double det1,det2;
      det1=inverse_mkeep( current_inv[alpha->lconf[i2]], m1, pars->NF[alpha->lconf[i2]] );
      cout << "(1) det_before: " << 1./det1;
*/
      r = update_row( &(adx[alpha->lconf[i1]][i2*pars->NF[alpha->lconf[i1]]]), row, current_inv[alpha->lconf[i1]], pars->NF[alpha->lconf[i1]] );
      current_x[alpha->lconf[i1]][row] = i2;

/*      if(abs(r)>SMALL) {
        det2=inverse_mkeep( current_inv[alpha->lconf[i2]], m, pars->NF[alpha->lconf[i2]] );
        cout << "; det_after: " << 1./det2;
        cout << "; r= " << r << "\n";
 
        if( (abs(det2)>SMALL && abs((det2*r - det1)/det2)>1e-5 ) ) {
          cout << "warn 1: d_matrix - old one " << alpha->lconf[i2] << " after swap of row " << row << " :\n";

          for(int i=0; i<pars->NF[alpha->lconf[i2]]; i++)
            for(int j=0; j<pars->NF[alpha->lconf[i2]]; j++)
              m[i][j] -= m1[i][j];
          write_m(m, pars->NF[alpha->lconf[i2]]);

          cout << "Replacement row:\n";
          for(int i=0; i<pars->NF[alpha->lconf[i2]]; i++) cout << adx[alpha->lconf[i2]][i2][i]-adx[alpha->lconf[i1]][i1][i] << ", ";
          cout << "\n";
        }
      }
*/
      if( abs(r)< SMALL ) {
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        return 0.;
      }
    }

    //find the row number of the second det
    for(row=0; (row<pars->NF[alpha->lconf[i2]])&&(current_x[alpha->lconf[i2]][row]!=i2); row++);
    
    if( row>=pars->NF[alpha->lconf[i2]] ) {
      cout << "Internal error i2 in swap\n";
      exit(-1);
    }

    if( ratio )
    {
      r *= det_ratio_row( &(adx[alpha->lconf[i2]][i1*pars->NF[alpha->lconf[i2]]]), row, current_inv[alpha->lconf[i2]], pars->NF[alpha->lconf[i2]] );

      //not swapping lconf here!
    }
    else
    {
/*      double det1,det2,r2;
      det1=inverse_mkeep( current_inv[alpha->lconf[i1]], m1, pars->NF[alpha->lconf[i1]] );
      cout << "(2) det_before: " << 1./det1;
*/
      double r2 = update_row( &(adx[alpha->lconf[i2]][i1*pars->NF[alpha->lconf[i2]]]), row, current_inv[alpha->lconf[i2]], pars->NF[alpha->lconf[i2]] );
//      r *= r2;
      current_x[alpha->lconf[i2]][row] = i1;

/*      if(abs(r2)>SMALL) {
        det2=inverse_mkeep( current_inv[alpha->lconf[i1]], m, pars->NF[alpha->lconf[i1]] );
        cout << "; det_after: " << 1./det2;
        cout << "; r= " << r2 << "\n";
        //inverse_mkeep(current_inv[alpha->lconf[i1]], m, pars->NF[alpha->lconf[i1]]);

        if( (abs(det1)>SMALL) && abs((det2*r2 - det1)/det1)>1e-5 ) {
          cout << "warn 2: d_matrix-old one " << alpha->lconf[i1] << " after swap of row " << row << " :\n";
          for(int i=0; i<pars->NF[alpha->lconf[i1]]; i++)
            for(int j=0; j<pars->NF[alpha->lconf[i1]]; j++)
              m[i][j] -= m1[i][j];
          write_m(m, pars->NF[alpha->lconf[i1]]);

          cout << "Replacement row:\n"; 
          for(int i=0; i<pars->NF[alpha->lconf[i1]]; i++) cout << adx[alpha->lconf[i1]][i1][i]-adx[alpha->lconf[i2]][i2][i] << ", ";
          cout << "\n";
        }
      }
*/
      wfn[alpha->lconf[i1]] *= r;
      wfn[alpha->lconf[i2]] *= r2;
      r *= r2;
      wf *= -r; //update the wf

      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
    }
    //destroy(m, pars->NF[0] );
  }


  //if( abs(r-1.) > .01 )
  //  cout << "swap (" << alpha->lconf[i1] << ", " << alpha->lconf[i2] << ") returning " << r << "\n";
  return -r;
}

//The crop function in not used in the U(1) state (flavor number conservation)
#if WFC
complex<double> u1::crop(int i1, int i2, bool ratio) {
#else
double u1::crop(int i1, int i2, bool ratio) {
#endif
  return .0;
}

void u1::find_starting_conf()
{
  int* Na0 = new int[NS];

  for(int n=0; n<NS; n++) {
    Na0[n] = N/NS;
  }

  //find_max_conf();
/*  int i, imax = 10000;
  for(i=0; i<imax; i++)
  {
    alpha->set_random_conf( Na0 );
    getwf();
    if( abs(wf)>SMALL ) break;
  }
  if( i==imax) cout << "ERROR: cannot find non-zero configuration...\n";
*/

  alpha->set_random_conf( Na0 );
  getwf();

  //cout << "Max config:\n";
  //alpha->print();
  //cout << "wf: " << wf << "\n";

  delete[] Na0;
}

void u1::find_max_conf() {}

int u1::insert_db()
{
  mysql_wrapper* wrapper = new mysql_wrapper();

  std::ostringstream os;

  os << "INSERT INTO paired (id, sites, lattice, txtdsc, NS, N0, ap, t1, t2, ";
  os << "P1, dP1, P2, dP2, P3, dP3, R4, dR4) VALUES (";
  os <<  "'', " << N << ", '" << alpha->mylattice->get_desc() << "', '" << pars->desc << "', " << NS << ", " << pars->NF[0] << ", " << (int)pars->ap;
  os << ", round(" << pars->t1 << ",3), round(" << pars->t2 << ",3), ";
  os << average[0] << ", " << sigma[0] << ", " << average[1] << ", " << sigma[1] << ", " << average[2] << ", " << sigma[2] << ", " << average[3] << ", " << sigma[3];
  os << ")";

  int res = wrapper->insert_qry( os.str() );

  delete wrapper;

  return res;
}

//corrects the wave fuction by pow(10,s 20), depending on the sign of s
void u1::correct_cff(bool s)
{
  double r;
  
  for(int n=0; n<NS; n++)
  {
    if( s )
      r = pow(10., +20./(double)pars->NF[n] );
    else
      r = pow(10., -20./(double)pars->NF[n] );
    normalize( r, n );
  }

  //cout << "TF: cff changed to " << cff[n_ext] << endl;
}

void u1::normalize(double r)
{
  for(int n=0; n<NS; n++) normalize(r, n);
}

void u1::normalize(double r, int n)
{
  if( n<0 || n>=NS ) {
    cout << "u1::normalize(): Invalid parameter...\n";
    exit(-1);
  }
  
  cout << "Normalizing " << n << " with r = " << r << endl;

  for(int j=0; j<N*pars->NF[n]; j++) adx[n][j] *= r;
  cff[n] *= r;
}
