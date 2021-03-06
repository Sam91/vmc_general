#include "u1hybrid.h"
#include "matrix.h"
#include <iomanip>

//define a general hybridized unpaired state on any lattice for NS flavors of fermions.
//argument l: linear system size
u1hybrid::u1hybrid(int l) : wavefunction( l )
{
  if( N%NS != 0 ) {
    cout << "N (number of sites " << N << ") must be a multiple of " << NS << "\n";
    exit(-1);
  }

#if WFC
  current_inv = new complex<double>[N*N];
  old_inv = new complex<double>[N*N];

  t = new (complex<double>(***[NS]));
  adx = new (complex<double>(*[NS]));
#if NP
  dadp = new (complex<double>(***[NP]));
  for(int np=0; np<NP; np++) dadp[np] = new (complex<double>(**[NS]));
#endif

  for(int n1=0; n1<NS; n1++)
  {
    adx[n1] =  new complex<double>[N*N];
    for(int p=0; p<NP; p++) dadp[p][n1] =  new complex<double>[N*N];

    t[n1] = new (complex<double>(**[NS]));
    for(int n2=0; n2<NS; n2++) {
      t[n1][n2] = createcomplex( N );
    }
  }
  d = createcomplex(NS, N);

#else
  current_inv = new double[N*N];
  old_inv = new double[N*N];

  t = new (double(***[NS]));
  adx = new (double(*[NS]));
#if NP
  dadp = new (double(**[NP]));
  for(int np=0; np<NP; np++) dadp[np] = new (double(*[NS]));
#endif

  for(int n1=0; n1<NS; n1++)
  {
    adx[n1] = new double[N*N];
    for(int p=0; p<NP; p++) dadp[p][n1] = new double[N*N];

    t[n1] = new (double(**[NS]));
    for(int n2=0; n2<NS; n2++) {
      t[n1][n2] = createdouble( N );
    }
  }
  d = createdouble(NS, N);

#endif

  for(int i=0; i<NS; i++) for(int j=0;j<N;j++) d[i][j]=0.;

  current_x = new int[N];
  old_x = new int[N];
  N0 = new double[NS];
  h = 0.;

  mu = new double[NS];

  e = new double[N*NS];
#if WFC
  v = createcomplex(N*NS);
#else
  v = createdouble(N*NS);
#endif
  Na0 = createdouble(NS*N, NS+2);
  occ = new int[N];
  for(int j=0; j<N; j++) occ[j]=0;
  isocc = new bool[NS*N];
  excit = -1;
}

u1hybrid::~u1hybrid()
{
  for(int n1=0; n1<NS; n1++)
  {
    for(int n2=0; n2<NS; n2++) destroy(t[n1][n2], N);
    delete[] t[n1];
    delete[] adx[n1];
    for(int p=0; p<NP; p++) {delete[] dadp[p][n1]; delete[] dadp[p];}
  }
  delete[] adx;
  delete[] t;
  destroy(d, NS);

#if NP
  delete[] dadp;
#endif

  delete[] current_inv;
  delete[] old_inv;

  delete[] current_x;
  delete[] old_x;
  delete[] N0;
  delete[] mu;

  delete[] e;
  destroy(v, N*NS);
  destroy(Na0, N*NS);
  delete[] occ;
  delete[] isocc;
}

/*void u1hybrid::print_t(int n1, int n2)
{
  int i1=0;
  for(int i=0; i<L; i++) {
    for(int j=0; j<L; j++) {
      cout << t[n1][n2][i1][] << ", ";
      i1++;
    }
    cout << "\n";
  }
}*/

void u1hybrid::set_excit(int x, int spin)
{
  cout << "Setting excit to "<< x << endl;
  if( x>= 0 ) {
    std::ostringstream os;
    os << bpars->desc << "-ex" << x << "_" << spin;
    bpars->desc = os.str();
  }
  excit = x;
  sp = spin;
}

//first, create the Hamiltonian matrix; We take (a x N + r) as indices
void u1hybrid::create_h0()
{
  for(int i1=0; i1<NS*N; i1++)
    for(int i2=0; i2<NS*N; i2++)
      v[i1][i2] = 0.;

  for(int n1=0; n1<NS; n1++) {
    for(int n2=0; n2<NS; n2++) {
      for(int i1=0; i1<N; i1++) {
        for(int i2=0; i2<N; i2++) {
          v[n1*N + i1][n2*N + i2] = t[n1][n2][i1][i2]; //hopping matrix
        }
#if WFC
        v[n1*N + i1][n2*N + i1] -= h*conj(d[n1][i1])*d[n2][i1];
#else
        v[n1*N + i1][n2*N + i1] -= h*d[n1][i1]*d[n2][i1];
#endif
      }
    }
  }
}

//calculate the MF spectrum and eigenvectors from the MF parameters
void u1hybrid::diagonalize()
{
  //construct the Hamiltonian matrix using t, d and, h
  create_h0(); //the Hamiltonian is saved in v

  //cout << "Hamiltonian matrix:\n";
  //write_m(v, N*NS);

  eigvects(v, e, N*NS); //e now contains the MF eigenvalues, and v the eigenvectors

/*  cout << "Spectrum: ";
  for(int i=0; i<N*NS; i++) cout << i << ":" << v[i] << "; ";
  cout << "\n";
*/

  //Go through all eigenstates and get their avg flavor number
  for(int ei=0; ei<NS*N; ei++) //loop through spectrum
  {
    Na0[ei][NS] = ei; //the last two values contain the spectral index and the energy
    Na0[ei][NS+1] = e[ei];
    for(int f=0; f<NS; f++) { //loop through flavors
      Na0[ei][f] = 0.;
      for(int j=0; j<N; j++) Na0[ei][f] += pow(abs(v[ei][f*N + j]), 2); //loop through sites
    }
  }
}

//determine the N out of the 2N eigenstates to be occupied in the GS; the indices of those states are stored in occ[]
void u1hybrid::construct_gs()
{
  double *Nstmp = new double[NS];

  if( abs(e[N]-e[N-1])<SMALL ) //degenerate at the Fermi level
  {
    cout << "WARN: Degenerate state; we need to fix this\n";
    cout << "(you may also consider changing boundary conditions or system size)\n";

    //get the lowest and highest index
    int ne0, ne1;
    for(ne0=0; ne0<N; ne0++) {
      if( abs(e[N-1] - e[ne0]) > SMALL ) occ[ne0] = ne0;
      else break;
    }
    for(ne1=ne0; ne1<N*NS; ne1++) {
      if( abs(e[N-1] - e[ne1]) > SMALL ) break;
    }
    cout << "Lowest and highest indices are " << ne0 << " and " << ne1 << ". ";
    cout << "We need to pick " << N-ne0 << " states out of " << ne1-ne0 << ".\n";

    double *Nstmp = new double{NS};
    for(int f=0; f<NS; f++) Nstmp[f] = (double)(N-ne0)/(double)NS;
    cout << "Nstmp: "; for(int f=0; f<NS; f++) cout << Nstmp[f] << ", "; cout << endl;

    for(int i=0; i<ne1-ne0; i++)
    {
      cout << i << ": " << e[i+ne0] << "; Na = ( ";
      for(int n=0; n<NS; n++) cout << Na0[i+ne0][n] << ", ";
      cout << ")\n";
    }
    cout << "Enter " << N-ne0 << " states: ";
    if( true ) //read in the states to use (choose by hand)
      for(int i=ne0; i<N; i++) {cin >> occ[i]; occ[i] += ne0;}

    else { //simply choose the first states that come
      int k = ne0;
      for(int i=ne0; i<ne1; i++)
      {
        bool cont = false;
        //cout << "Ntmp: "; for(int f=0; f<NS; f++) cout << Nstmp[f] << ", "; cout << endl;
        for(int f=0; f<NS; f++) {if( Nstmp[f]-Na0[i][f] < -.1 ) cont = true;}
        if( cont ) continue;

        for(int f=0; f<NS; f++) Nstmp[f] -= Na0[i][f];

        cout << i-ne0 << " ";
        if( k>= N) break;
        occ[k] = i; k++;
      }
      cout << "\n";
    }
    delete[] Nstmp;
  } else { //non-degenerate case
    for(int i=0; i<N; i++) occ[i] = i;
  }

  //get the occupied states
  for(int ie=0; ie<NS*N; ie++) isocc[ie] = false;
  for(int j=0; j<N; j++) isocc[ occ[j] ] = true;

  delete[] Nstmp;
}

//Construct a spinfull excitation on top of the GS

//full bandwidth
void u1hybrid::construct_ex3()
{
  //first, find the index if the lowest down spin
  int ie;
  for(ie=0; ie<N; ie++) {
    cout << "e[" << ie << "] = " << e[ie] << ": Na[0] = " << Na0[ie][0] << endl;
    if( Na0[ie][0]>.8 ) break;
  }

  //find highest up excit
  int iex;
  for(iex=N*NS-1; iex>=0; iex--) {
    if( isocc[iex] ) continue;
    cout << "e[" << iex << "] = "<< e[iex] << ": Na[0] = " << Na0[iex][0] << endl;
    if( Na0[iex][sp]>.8 ) break;
  }

  //replace the state
  for(int i=0; i<N; i++) 
    if( occ[i]==ie ) {occ[i] = iex; break;}
  isocc[ie] = false; isocc[iex] = true;
}

//half bandwidth
void u1hybrid::construct_ex1()
{
  //find highest down spin hole
  int ie;
  for(ie=NS*N-1; ie>=0; ie--) {
    if( !isocc[ie] ) continue;
    cout << "e[" << ie << "] = " << e[ie] << ": Na[0] = " << Na0[ie][0] << endl;
    if( Na0[ie][0]>.8 ) break;
  }

  //find highest up spin excit
  int iex;
  for(iex=N*NS-1; iex>=0; iex--) {
    if( isocc[iex] ) continue;
    cout << "e[" << iex << "] = "<< e[iex] << ": Na[0] = " << Na0[iex][0] << endl;
    if( Na0[iex][sp]>.8 ) break;
  }

  //replace the state
  for(int i=0; i<N; i++)
    if( occ[i]==ie ) {occ[i] = iex; break;}
  isocc[ie] = false; isocc[iex] = true;
}

void u1hybrid::construct_ex2()
{
  //find the lowest down spin hole
  int ie;
  for(ie=0; ie<NS*N; ie++) {
    if( !isocc[ie] ) continue;
    cout << "e[" << ie << "] = " << e[ie] << ": Na[0] = " << Na0[ie][0] << endl;
    if( Na0[ie][0]>.8 ) break;
  }

  //find lowest up spin excit
  int iex;
  for(iex=0; iex<NS*N; iex++) {
    if( isocc[iex] ) continue;
    cout << "e[" << iex << "] = "<< e[iex] << ": Na[0] = " << Na0[iex][0] << endl;
    if( Na0[iex][sp]>.8 ) break;
  }

  //replace the state
  for(int i=0; i<N; i++)
    if( occ[i]==ie ) {occ[i] = iex; break;}
  isocc[ie] = false; isocc[iex] = true;
}

//Fermi surface excitations
void u1hybrid::construct_ex0()
{
  //find the highest down spin hole
  int ie;
  for(ie=NS*N-1; ie>=0; ie--) {
    if( !isocc[ie] ) continue;
    cout << "e[" << ie << "] = " << e[ie] << ": Na[0] = " << Na0[ie][0] << endl;
    if( Na0[ie][0]>.8 ) break;
  }

  //find lowest up spin excit
  int iex;
  for(iex=0; iex<NS*N; iex++) {
    if( isocc[iex] ) continue;
    cout << "e[" << iex << "] = "<< e[iex] << ": Na[0] = " << Na0[iex][0] << endl;
    if( Na0[iex][sp]>.8 ) break;
  }

  //replace the state
  for(int i=0; i<N; i++)
    if( occ[i]==ie ) {occ[i] = iex; break;}
  isocc[ie] = false; isocc[iex] = true;
}

//initialize the hybridized U(1) state (ground state)
void u1hybrid::create_ad()
{
  cout << "Entering u1hybrid::create_ad()\n";

  //construct Hamiltonian and diagonalize it
  diagonalize();

  //construct the MF GS
  construct_gs();
  cout << "occ: "; for(int i=0; i<N; i++) cout << occ[i] << ", "; cout << endl;

  //construct an exitation on top of the gs

  if( excit >= 0 )
  {
    if( excit==0 ) construct_ex0();
    else if(excit==1) construct_ex1();
    else if(excit==2) construct_ex2();
    else if(excit==3) construct_ex3();
    cout << "occ: "; for(int i=0; i<N; i++) cout << occ[i] << ", "; cout << endl;
  }
  //cout << "occ: "; for(int i=0; i<N; i++) cout << occ[i] << ", "; cout << endl;

  //construct an exitation on top of the gs
  //construct_ex();
  //cout << "occ: "; for(int i=0; i<N; i++) cout << occ[i] << ", "; cout << endl;

/*
    //calculate the A-matrix for the derivative
    for(int i1=0; i1<N*NS; i1++)
      for(int i2=0; i2<N*NS; i2++)
        for(int k1=0; k1<N*NS; k1++)
          for(int k2=0; k2<N*NS; k2++) {
            dad[i1][i2][k1][k2] = 0.;
            for(int n1=0; n1<N; n1++) {
              for(int n2=N; n2<N*NS; n2++) {
#if WFC
                dad[i1][i2][k1][k2] += conj(v[n1][i1])*v[n1][i2]*v[n2][k1]*conj(v[n2][k2])/(e[n1] - e[n2]);
#else
                dad[i1][i2][k1][k2] += v[n1][i1]*v[n1][i2]*v[n2][k1]*v[n2][k2]/(e[n1] - e[n2]);
#endif
              }
            }
          }
*/

  //write the chosen N eigenvectors into the adx matrix and calculate the avg flavor number
  for(int f=0; f<NS; f++) N0[f] = 0.;

  for(int ie=0; ie<N; ie++)
  {
    for(int j=0; j<N; j++) {
      for(int f=0; f<NS; f++)
      {
        adx[f][j*N+ie] = v[ occ[ie] ][f*N + j]; //note that the eigenvalue index is the last index in adx!
        N0 [f] += pow( abs(v[ occ[ie] ][f*N + j]), 2);
      }
    }
  }

  cout << "Avg flavor number before projection: ";
  for(int f=0; f<NS; f++) cout << N0[f] << ", "; cout << "\n";

  //set NF (wf projection) to the mean-field value
  for(int f=0; f<NS; f++) NF[f] = round(N0[f]);

  //cout << "abs diff " << abs(v[i]-v[i-1]) << "\n";
/*  if( abs(v[ne]-v[ne-1]) < 1e-8 )
    cout << "WARN: degenerate state\n";
  cout << "Fermi energy is " << v[ne-1] << ", " << v[ne] << " at position " << ne-1 << "\n"; 
*/

  cff = 1.;
  normalize( 3. );

/*  for(int n=0; n<NS; n++)
  {
    cout << "adx[ " << n << " ]\n";
    write_m(adx[n],N);
  }*/
}

void u1hybrid::set_h(double h0) {cout << "Setting h to " << h0 << "\n"; h = h0;}

//this function tries to find an apropriate chemical potential mu[0] such that the avg flavor number before projections has N[0]=N[1]
int u1hybrid::findmu()
{
  int n, nmax;
  cout << "Entering findmu();\n";

  double dm = .1; //delta mu

  nmax = 10;
  for(n=0; n<10; n++)
  {
    create_ad();

    if( abs(N0[0]-N0[1])<.1 ) {
      cout << "Appropriate muz is " << mu[0] << "\n";
      break;
    } else {
      if( N0[0]>N0[1] ) mu[0] -= dm;
      else mu[0] += dm;
    }
    for(int i=0; i<N; i++)
      for(int n=0; n<NS; n++)
        t[n][n][i][i] = -mu[n];
  }

  if(n==nmax) {
    cout << "No apropriate mu_z found...\n";
    return -1;
  }
  return 0;
}

//set a 0-flux t-t' hopping state on some lattice torus (make sure the lattice has been set)
//set a hopping that is diagonal in flavor index
void u1hybrid::set_hopping(double tt1, double tt2, bool ap)
{
  for(int n=0; n<NS; n++) mu[n] = 0.;

  int i2;

  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n1=0; n1<NS; n1++)
        for(int n2=0; n2<NS; n2++) t[n1][n2][i][j] = 0.;

  for(int i1=0; i1<N; i1++)
  {
    for(int n=0; n<NS; n++)
    {
      for(int k=1; k<=alpha->mylattice->connectivity[i1][0][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][0][k];
        if( ap && ( abs(i1%L - i2%L)>L/2 ) ) t[n][n][i1][i2] += -tt1;
        else t[n][n][i1][i2] += tt1;
      }
      for(int k=1; k<=alpha->mylattice->connectivity[i1][1][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][1][k];
        if( ap && ( abs(i1%L - i2%L)>L/2 ) ) t[n][n][i1][i2] += tt2;
        else t[n][n][i1][i2] += -tt2;
      }
      //add the chemical potential
//      t[n][n][i1][i1] -= mu[n];
    }
  }
}

//set a special hopping pattern on the triangular lattice
void u1hybrid::set_hopping3(double* tt1, double* tt2, bool ap)
{
  for(int n=0; n<NS; n++) mu[n] = 0.;

  int i2;

  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int n1=0; n1<NS; n1++)
        for(int n2=0; n2<NS; n2++) t[n1][n2][i][j] = 0.;

  for(int i1=0; i1<N; i1++)
  {
    for(int n=0; n<NS; n++)
    {
      for(int k=1; k<=alpha->mylattice->connectivity[i1][0][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][0][k];
        if( ap && ( abs(i1%L - i2%L)>L/2 ) ) t[n][n][i1][i2] += -tt1[k-1];
        else t[n][n][i1][i2] += tt1[k-1];
      }
      for(int k=1; k<=alpha->mylattice->connectivity[i1][1][0]; k++) {
        i2 = alpha->mylattice->connectivity[i1][1][k];
        if( ap && ( abs(i1%L - i2%L)>L/2 ) ) t[n][n][i1][i2] += tt2[k-1];
        else t[n][n][i1][i2] += -tt2[k-1];
      }
    }
  }
}

int u1hybrid::findmu( int nz ) //tries to find a mu[0] such that Nz = nz before projection
{
  int n, nmax;
  cout << "Entering findmu(" << nz << ");\n";

  NF[0] = nz;
  for(int i=1; i<NS; i++) NF[i] = (N-nz)/(NS-1);

  double dm = .05; //delta mu

  double delta_n, last_mu2, last_mu, last_d2, last_d;
  last_mu = mu[0]-1.;

  nmax = 30;
  for(n=0; n<nmax; n++)
  {
    cout << " muz = " << mu[0] << "\n";
    create_ad();

    delta_n = abs(N0[0]-nz);
    if( delta_n<.1 )
    {
      cout << "Appropriate muz is " << mu[0] << "\n";
      break;
    } else
    {
      last_mu2 = last_mu; last_mu = mu[0]; //save the last muz
      last_d2 = last_mu2; last_d = delta_n;

      if( N0[0]>nz ) mu[0] -= dm;
      else mu[0] += dm;

      if( abs(last_mu2-mu[0])<1e-10 )
      {
        cout << "Switching between two mus: ";
        if(last_d2 < last_d)
          mu[0] = last_mu2;
        else
          mu[0] = last_mu;
        cout << "Let's take muz = " << mu[0] << "\n";
        break;
      }
    }
    for(int i=0; i<N; i++)
      for(int n=0; n<NS; n++)
        t[n][n][i][i] = -mu[n];
  }

  if(n==nmax) {
    cout << "No apropriate muz found...\n";
    return -1;
  }
  return 0;
}

void u1hybrid::backup_data()
{
  memcpy( old_x, current_x, N*sizeof(int) );

#if WFC
  memcpy( old_inv, current_inv, N*N*sizeof(complex<double>) );
#else
  memcpy( old_inv, current_inv, N*N*sizeof(double) );
#endif

  wf_old = wf;
}

void u1hybrid::restore_data()
{
  int* tmp;
#if WFC
  complex<double>* temp;
#else
  double* temp;
#endif

  tmp = current_x; current_x = old_x; old_x = tmp;
  temp = current_inv; current_inv = old_inv; old_inv = temp;

  wf = wf_old;
}

//get the wf = <alpha|psi>
void u1hybrid::getwf()
{  
  //loop over all sites and pick the adx eigenvector component corresponding to the flavor on that site
  for(int j=0; j<N; j++)
  {
#if WFC
    //memcpy(d_matrix[j], &(adx[ alpha->lconf[j] ][j]), N*sizeof(complex<double>) );
    memcpy(&(current_inv[j*N]), &(adx[ alpha->lconf[j] ][j*N]), N*sizeof(complex<double>) );
#else
    memcpy(&(current_inv[j*N]), &(adx[ alpha->lconf[j] ][j*N]), N*sizeof(double) );
#endif
    //for(int n=0; n<N; n++) d_matrix[j*N+n] = adx[ alpha->lconf[j] ][j*N+n];
    current_x[j] = j;
  }

  //cout << "d_matrix:\n";
  //write_m(d_matrix, N);

  //wf = inverse_m(d_matrix, current_inv, N );
  wf = inverse(current_inv, N );

  if( abs(wf)>1e300 ) {
    cout << "Error: divergent wf\n";
    exit(-1);
  }
  //cout << "wf norm: " << wf << "\n";

#if JASTROW
  wf *= jastrow3();
#endif

  //cout << "getwf(): " << std::scientific << wf << "\n";
}

//calculate the gradient of the wave function with respect to the variational parameters
//(note that the inverse matrix for the current configuration should be updated; the differential matrix is calculated from scratch)
void u1hybrid::get_dwf()
{
#if NP
  //loop over variational parameters
  for(int p=0; p<NP; p++)
  {
    dwf[p] = 0.;
    //loop over all sites and pick the dadp eigenvector component corresponding to the flavor on that site
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++)
        dwf[p] += current_inv[k*N+j] * dadp[p][ alpha->lconf[j] ][j*N+k];
  }

  for(int i1=0; i1<N*NS; i1++)
    for(int i2=0; i2<N*NS; i2++)
      dF[i1][i2] = 0.;

  for(int i1=0; i1<N; i1++)
  {
    n1 = alpha->lconf[i1];
    k1=n1*N+i1;

    for(int i2=0; i2<N; i2++) {
      if(i1==i2 && n1==n2) {
        dF[n1*N + i1][n1*N + i1] = wf;
        continue;
      }
      for(int n2=0; n2<NS; n2++) {
        k2=n2*N+i2;
        for(int n=0; n<N; n++) dF[k1][k2] += current_inv[i1*N+n]*adx[n][n2*N+i2];
      }
    }
  }
#endif
}

//swap operator of two sites i1 and i2 (corresponding to virtual_replacement)
//If ratio=True, we only return the ratio of the wf after the swap, without updating the current wf
#if WFC
complex<double> u1hybrid::swap(int i1, int i2, bool ratio)
{
  complex<double> r;
#else
double u1hybrid::swap(int i1, int i2, bool ratio)
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
    if( abs(wf)<SMALL )
    {
      if( ratio ) { cout << "ERROR: incorrect ratio for vanishing wf! (wf="<< wf << ")\n"; exit(-1);}

      cout << "vanishing wf prior to swap... calculating from scratch ["<< wf << "]\n";

      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
      wf_old = wf;
      getwf(); //uses the current alpha configuration (needs to be swapped at that point)
      return wf/wf_old; //swap will be accepted
    }

/*    cout << "swap(" << i1 << ", " << i2 << ", " << ratio << ")\n";
    cout << "alpha: ";
    for(int i=0; i<N; i++) cout << alpha->lconf[i] << ", ";
    cout << "\n";
*/
/*    double **m = createdouble( N ); double **m1 = createdouble( N );
    double det1, det2;
    det1=inverse_mkeep( current_inv, m1, N );
*/
/*    //cout << "(1) det_before: " << 1./det1;

    //cout << "Current d:\n";
    //write_m(m1, N);
*/
/*    cout << "We are going to replace row " << i1 << " by {";
    for(int i=0; i<N; i++) cout << adx[ alpha->lconf[i2] ][i1][i] << ", ";
    cout << "}\n";
    cout << " and row " << i2 << " by {";
    for(int i=0; i<N; i++) cout << adx[ alpha->lconf[i1] ][i2][i] << ", ";
    cout << "}\n";
*/
    int row1;
    for(row1=0; (row1<N)&&(current_x[row1]!=i1); row1++);
    if( row1>=N ) {
      cout << "Internal error i1 in swap\n";
      exit(-1);
    }

    if( ratio )
    {
      backup_data(); //For hybridized wf, we need to backup
      r = update_row( &(adx[ alpha->lconf[i1] ][i2*N]), row1, current_inv, N );

      if( abs(r)<SMALL )
      { //For hybridized states, we cannot return here. We need to calculate from scratch
//        cout << "(p) Intermediate zero; recalculating from scratch\n";
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        getwf();
#if WFC
        complex<double> rr = wf/wf_old;
#else
        double rr = wf/wf_old;
#endif
        restore_data();
        c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        return rr;
      } //else cout << "(p) Intermediate non-zero\n";

    } else {

      //cout << "Current_inv:\n";
      //write_m(current_inv, N);

      r = update_row( &(adx[ alpha->lconf[i1] ][i2*N]), row1, current_inv, N ); //this updates the current_inv matrix
/*
      if( abs(r)>SMALL ) {
        det2=inverse_mkeep( current_inv, m, N );
        //cout << "; det_after: " << 1./det2;
        //cout << "; r= " << r << "\n";
 
        if( (abs(det2)>SMALL && abs((det2*r - det1)/det2)>1e-5 ) )
        {
          cout << "warn 1: d_matrix - old one " << alpha->lconf[i2] << " after swap of row " << i1  << "(det1=" << det1 << "; det2=" << det2 << "; r=" << r << "):\n";

          for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
              m[i][j] -= m1[i][j];
          write_m(m, N);

          cout << "Replacement row:\n";
          for(int i=0; i<N; i++) cout << adx[alpha->lconf[i2]][i2*N+i]-adx[alpha->lconf[i1]][i1*N+i] << ", ";
          cout << "\n";
        }
      } else cout << "\n";
*/
      if( abs(r)<SMALL )
      {
        int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
        //cout << "Intermediate zero; recalculating from scratch\n";
        getwf();
        
        //destroy(m, N); destroy(m1, N);
        if( abs(wf_old)>SMALL )
          return wf/wf_old;
        else return 1.;
      } //else cout << "Intermediate non-zero\n";
    }

    int row2;
    for(row2=0; (row2<N)&&(current_x[row2]!=i2); row2++);
    if( row2>=N ) {
      cout << "Internal error i2 in swap\n";
      exit(-1);
    }

#if JASTROW
    r *= jastrow3(i1, i2);
#endif

    if( ratio )
    {
      r *= det_ratio_row( &(adx[ alpha->lconf[i2] ][i1*N]), row2, current_inv, N );

      restore_data();
      //not swapping lconf here!
    }
    else
    {
/*      double det1,det2,r2;
      det1=inverse_mkeep( current_inv, m1, N );
      cout << "(2) det_before: " << 1./det1;
*/
      r *= update_row( &(adx[ alpha->lconf[i2] ][i1*N]), row2, current_inv, N );
      //r *= r2;

      current_x[row1] = i2;
      current_x[row2] = i1;

/*      if(abs(r2)>SMALL) {
        det2=inverse_mkeep( current_inv, m, N );
        //cout << "; det_after: " << 1./det2;
        //cout << "; r= " << r2 << "\n";

        if( (abs(det1)>SMALL) && abs((det2*r2 - det1)/det1)>1e-5 )
        {
          cout << "warn 2: d_matrix-old one " << alpha->lconf[i1] << " after swap of row " << i2 << "(det1=" << det1 << "; det2=" << det2 << "; r=" << r2 << "):\n";
          for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
              m[i][j] -= m1[i][j];
          write_m(m, N);

          cout << "Replacement row:\n"; 
          for(int i=0; i<N; i++) cout << adx[alpha->lconf[i1]][i1][i]-adx[alpha->lconf[i2]][i2][i] << ", ";
          cout << "\n";
        }
      } else cout << "\n";
*/
      wf *= -r; //update the wf

      //swap in lconf
      int c = alpha->lconf[i1]; alpha->lconf[i1] = alpha->lconf[i2]; alpha->lconf[i2] = c;
    }
//    destroy(m, N ); destroy(m1, N);
  }

  //if( abs(r-1.) > .01 )
  //  cout << "swap (" << alpha->lconf[i1] << ", " << alpha->lconf[i2] << ") returning " << r << "\n";
  return -r;
}

//The crop function in not used in the U(1) state (flavor number conservation)
#if WFC
complex<double> u1hybrid::crop(int i1, int i2, bool ratio) {
#else
double u1hybrid::crop(int i1, int i2, bool ratio) {
#endif
  return .0;
}

void u1hybrid::find_max_conf() {}

//correct the wf normalization by pow(10,s 20), depending on the sign of s
void u1hybrid::correct_cff(bool s)
{
  double r;
  
  for(int n=0; n<NS; n++)
  {
    if( s )
      r = pow(10., +40./(double)N );
    else
      r = pow(10., -60./(double)N );
  }

  normalize( r );
  //cout << "TF: cff changed to " << cff[n_ext] << endl;
}

void u1hybrid::normalize(double r)
{
  cout << "Normalizing with r = " << r << endl;

  for(int n=0; n<NS; n++)
    for(int j=0; j<N*N; j++) adx[n][j] *= r;
  cff *= r;
}

