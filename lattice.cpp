#include "lattice.h"
#include "helperfunctions.h"

lattice::lattice( int L0 )
{
  this->L = L0;

  this->Q = SUBL;

  LD = pow(L,DIM); N = LD*SUBL;

  cout << "Creating lattice with L=" << L << " and " << "Q=" << Q << "." << endl;

  //max number of neighbors we want to store in the links/connectivity matrices
  nbr=10;

  /* The site relations are defined in the following way
   * plaquette[site index j1][type of plaquette][number of plaquettes on this site, j2, j3, j4, ...]
  */

  connectivity = new (int(**[ this->N ])); //lets take first and second neighbor
  links = new (int(**[ this->N ]));
  plaquette4 = new (int(**[ this->N ])); //four-site plaquettes

  for(int i=0; i<N; i++)
  {
    connectivity[i] = createint(nbr, 10);
    links[i] = createint(nbr, 10);
    plaquette4[i] = createint(nbr, 20);
  }
  n1 = new int[DIM];
  n2 = new int[DIM];
}

lattice::~lattice()
{
  for(int i=0; i<N; i++)
  {
    destroy(connectivity[i], nbr);
    destroy(links[i], nbr);
    destroy(plaquette4[i], nbr);
  }
  delete[] connectivity;
  delete[] links;
  delete[] plaquette4;

  delete[] n1; delete[] n2;
}

std::string lattice::get_desc() {return desc;}

// mapping from cartesian coordinates to linear coordinates
int lattice::j(int* n, int q)
{
  //move to the periodic torus
  for(int d=0; d<DIM; d++)
  {
    if( n[d]>= L ) n[d] -= L;
    if( n[d] < 0 ) n[d] += L;
  }

#if DIM == 1
#if SUBL>1
  return SUBL*n[0] + q;
#else
  return n[0];
#endif

#elif DIM == 2
#if SUBL>1
    return SUBL*(n[0] + n[1]*L) + q;
#else
    return n[0] + n[1]*L;
#endif

#else
  cout << "ERROR: Undefined dimension\n";
  exit(-1);

#endif
}

// inverse mapping, from linear to cartesian coordinates
void lattice::getnq(int* n, int &q, int jj)
{
  //cout << "L=" << L << endl;
  //cout << "this->L=" << this->L << endl;

  if(jj<0 || jj>=N) {
    cout << "ERROR: undefined mod...\n";
    exit(-1);
  }

#if SUBL>1
  q = jj % SUBL;
  int j0 = jj / SUBL;

#else
  int j0 = jj;

#endif

#if DIM == 1
  n[0] = j0;

#elif DIM == 2
  n[0] = j0 % L;
  n[1] = j0 / L;

#else
  cout << "ERROR: undefined dimension\n";
  exit(-1);

#endif
}

//restrain the values to the LxL using periodc boundary conditions
void lattice::torus(int* n)
{
  for(int d=0; d<DIM; d++)
  {
    if( n[d]>= L ) n[d] -= L;
    if( n[d] < 0 ) n[d] += L;
  }
}

//sets the connectivity and links matrices to a kagome (currently only nearest neighbor!)
void lattice::set_kagome()
{
  if( DIM != 2 )
  {
    cout << "lattice::set_kagome(): need to compile with DIM = 2." << endl;
    exit(-1);
  }
  if( SUBL != 3 )
  {
    cout << "lattice::set_kagome(): need to compile with SUBL = 3." << endl;
  }

  int n;
  int nx, ny, q;
  int q2, i2;

  //checking for consistency
  desc = "kagome";
  cout << "Setting lattice to " << desc << ".\n";

/*  if( Q != 3 )
  {
    cout << "lattice::set_kagome(): cannot set kagome; number of sites in the unit cell Q is different from 3 "<<endl;
    exit(-1);
  }
*/
/*
  //loop over all UNIT CELLS of the lattice (cell: three sites of the hexagon) 
  for(int nx=0; nx<L; nx++) // x-coordinate
  {
    for(int ny=0; ny<L; ny++) // y-coordinate
    {
      for(int q=0; q<Q; q++) // loop over positions in this unit cell
      {
*/

  //loop over all sites of the lattice
  for(int i1=0; i1<N; i1++)
  {
    getnq(n1, q, i1); nx = n1[0]; ny = n1[1];

//        n1[0] = nx; n1[1] = ny;
//        i1 = j(n1, q);

        //setting nearest neighbors: links =  distinct links in the direction 0, 60 and 120 degree; connectivity = all links
        n = 0;
        links[i1][n][0] = 2; //number of nn for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        { //ordering: straight lignes, then one rotation to the right
          // +y
          n2[0] = nx; n2[1] = ny; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // -y
          n2[0] = nx; n2[1] = ny-1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +x+y
          n2[0] = nx+1; n2[1] = ny; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;
          
          // -x-y
          n2[0] = nx; n2[1] = ny-1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }  

        if( q==1 ) // second site in the unit cell
        {
          // -x
          n2[0] = nx; n2[1] = ny; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x
          n2[0] = nx+1; n2[1] = ny; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +y
          n2[0] = nx; n2[1] = ny+1; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;

          // -y
          n2[0] = nx; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==2 ) // third site in the unit cell
        {
          // -x-y
          n2[0] = nx-1; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x+y
          n2[0] = nx; n2[1] = ny+1; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -x
          n2[0] = nx-1; n2[1] = ny; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;

          // +x
          n2[0] = nx; n2[1] = ny; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        // second neighbors
        n = 1;
        links[i1][n][0] = 2;//number of 2nd for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        { //ordering: straight lines, then one rotation to the left

          // -x+y
          n2[0] = nx; n2[1] = ny; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x-y
          n2[0] = nx+1; n2[1] = ny-1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -2x-y
          n2[0] = nx-1; n2[1] = ny-1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;

          // +2x+y
          n2[0] = nx+1; n2[1] = ny; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==1 ) // second site in the unit cell
        {
          // -2x-y
          n2[0] = nx-1; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x+y
          n2[0] = nx+1; n2[1] = ny+1; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -x-2y
          n2[0] = nx; n2[1] = ny-1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;

          // +x+2y
          n2[0] = nx+1; n2[1] = ny+1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==2 ) // third site in the unit cell
        {
          // -x-2y
          n2[0] = nx-1; n2[1] = ny-1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x+2y
          n2[0] = nx; n2[1] = ny+1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +x-y
          n2[0] = nx; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;

          // -x+y
          n2[0] = nx-1; n2[1] = ny+1; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

/* deactivating non-diagonal third neighbor for now, because we are not using it for now
        // third neighbors (non-diagonal)
        n = 2;
        links[i1][n][0] = 2;  //number of 2nd for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        {
          // +x,0
          n2[0] = nx+1; n2[1] = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +x+y,0
          n2[0] = nx+1; n2[1] = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x,0
          n2[0] = nx-1; n2[1] = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;
          // -x-y,0
          n2[0] = nx-1; n2[1] = ny-1; q2 = 0; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==1 ) // second site in the unit cell
        {
          // +x,1
          n2[0] = nx+1; n2[1] = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +y,1
          n2[0] = nx; n2[1] = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x,1
          n2[0] = nx-1; n2[1] = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;
          // -y,1
          n2[0] = nx; n2[1] = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==2 ) // third site in the unit cell
        {
          // +x+y,2
          n2[0] = nx+1; n2[1] = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +y,2
          n2[0] = nx; n2[1] = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x-y,2
          n2[0] = nx-1; n2[1] = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][3] = i2;
          // -y,2
          n2[0] = nx; n2[1] = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][4] = i2;
        }
*/
        // diagonal neighbors (inside the hexagon)
        n = 2;
        links[i1][n][0] = 1;//number of 2nd for a given site
        connectivity[i1][n][0] = 2;

        if( q==0 ) // first site in the unit cell
        {
          // -2x
          n2[0] = nx-1; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x
          n2[0] = nx+1; n2[1] = ny; q2 = 0; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][2] = i2;
        }

        if( q==1 ) // second site in the unit cell
        {
          // -2x-2y
          n2[0] = nx-1; n2[1] = ny-1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x+2y
          n2[0] = nx+1; n2[1] = ny+1; q2 = 1; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][2] = i2;
        }

        if( q==2 ) // third site in the unit cell
        {
          // -2y
          n2[0] = nx; n2[1] = ny-1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2y
          n2[0] = nx; n2[1] = ny+1; q2 = 2; //torus(nx2, ny2);
          i2 = j(n2, q2); connectivity[i1][n][2] = i2;
        }

//      } //loop over sites in unit cell
//    } //loop over sites (x)
//  } //loop over sites (y)

  } //loop over all sites i1
}

//sets the connectivity and links matrices to a triangular lattice (currently up to third neighbor, and four-site plaquette term)
void lattice::set_triangular()
{
  if( DIM != 2 )
  {
    cout << "lattice::set_triangular(): need to compile with DIM=2." << endl;
    exit(-1);
  }
  if( SUBL != 1)
  {
    cout<< "lattice::set_triangular(): need to compile with SUBL=1." << endl;
    exit(-1);
  }

  int n;
  int i2;
  int nx, ny, q;

  desc = "triangular";
  cout << "Setting lattice to " << desc << ".\n";

  for(int i1=0; i1<N; i1++)
  {
    getnq(n1, q, i1); nx = n1[0]; ny = n1[1];

    //nearest neighbor
    n=0;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //x
    n2[0] = nx+1; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //x+y
    n2[0] = nx+1; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //y
    n2[0] = nx; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-x
    n2[0] = nx-1; n2[1] = ny; //torus(nx2, ny2);
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //-x-y
    //y direction
    n2[0] = nx-1; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][5] = i2;

    //-y
    n2[0] = nx; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][6] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //2x+y
    n2[0] = nx+2; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //x+2y
    n2[0] = nx+1; n2[1] = ny+2; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x+y
    n2[0] = nx-1; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-2x-y
    n2[0] = nx-2; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //-x-2y
    n2[0] = nx-1; n2[1] = ny-2; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][5] = i2;

    //x-y
    n2[0] = nx+1; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][6] = i2;

    //3rd neighbor
    n=2;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //2x
    n2[0] = nx+2; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //2x+2y
    n2[0] = nx+2; n2[1] = ny+2; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //2y
    n2[0] = nx; n2[1] = ny+2; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-2x
    n2[0] = nx-2; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //-2x-2y
    n2[0] = nx-2; n2[1] = ny-2; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][5] = i2;

    //-2y
    n2[0] = nx; n2[1] = ny-2; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][6] = i2;

    //4-site plaquette term
    n=0;
    plaquette4[i1][n][0] = 3;
    n2[0] = nx+1; n2[1] = ny  ; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][1] = i2;
    n2[0] = nx+1; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][2] = i2;
    n2[0] = nx  ; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][3] = i2;

    n2[0] = nx+1; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][4] = i2;
    n2[0] = nx+1; n2[1] = ny+2; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][5] = i2;
    n2[0] = nx  ; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][6] = i2;

    n2[0] = nx+1; n2[1] = ny  ; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][7] = i2;
    n2[0] = nx+2; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][8] = i2;
    n2[0] = nx+1; n2[1] = ny+1; /*torus(nx2, ny2);*/ i2 = j(n2, 0); plaquette4[i1][n][9] = i2;
  }
}

//square lattice on a torus
void lattice::set_square()
{
  if( DIM != 2 )
  {
    cout << "lattice::set_square(): need to compile with DIM = 2." << endl;
    exit(-1);
  }
  if( SUBL != 1)
  {
    cout << "lattice::set_square(): need to compile with SUBL = 1." << endl;
    exit(-1);
  }

  int n, i2;
  int nx, ny, q;

  desc = "square";
  cout << "Setting lattice to " << desc << ".\n";

  for(int i1=0; i1<LD; i1++)
  {
    getnq(n1, q, i1); nx = n1[0]; ny = n1[1];

    //nearest neighbor
    n=0;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x
    n2[0] = nx+1; n2[1] = ny; //torus(nx2, ny2);
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //y
    n2[0] = nx; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x
    n2[0] = nx-1; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][3] = i2;

    //-y
    n2[0] = nx; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x+y
    n2[0] = nx+1; n2[1] = ny+1; //torus(nx2, ny2);
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //-x+y
    n2[0] = nx-1; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x-y
    n2[0] = nx-1; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][3] = i2;

    //x-y
    n2[0] = nx+1; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //4-site plaquette term
    n=0;
    plaquette4[i1][n][0] = 1;
    n2[0] = nx+1; n2[1] = ny  ; /*torus(nx2, ny2)*/; i2 = j(n2, 0); plaquette4[i1][n][1] = i2;
    n2[0] = nx+1; n2[1] = ny+1; /*torus(nx2, ny2)*/; i2 = j(n2, 0); plaquette4[i1][n][2] = i2;
    n2[0] = nx  ; n2[1] = ny+1; /*torus(nx2, ny2)*/; i2 = j(n2, 0); plaquette4[i1][n][3] = i2;
  }
}

void lattice::set_checkerboard()
{
  if( DIM != 2 )
  {
    cout << "lattice::set_checkerboard(): need to compile with DIM = 2." << endl;
    exit(-1);
  }
  if( SUBL != 2 )
  {
    cout << "lattice::set_checkerboard(): need to compile with SUBL = 2." << endl;
    exit(-1);
  }

  int n, i2;

  if(L%2!=0) {
    cout << "Error: cannot set checkerboard... (L=" << L << ")\n";
    exit(-1);
  }

  desc = "checkerboard";
  cout << "Setting lattice to " << desc << ".\n";

  int nx, ny, q;

  for(int i1=0; i1<LD; i1++)
  {
    getnq(n1, q, i1); nx = n1[0]; ny = n1[1];

    //nearest neighbor
    n=0;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x
    n2[0] = nx+1; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

    //y
    n2[0] = nx; n2[1] = ny+1; //torus(nx2, ny2); 
    i2 = j(n2, 0); links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x
    n2[0] = nx-1; n2[1] = ny; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][3] = i2;

    //-y
    n2[0] = nx; n2[1] = ny-1; //torus(nx2, ny2); 
    i2 = j(n2, 0); connectivity[i1][n][4] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 1;
    connectivity[i1][n][0] = 2;

    if( ((nx+ny)%2)==0 )
    { //right

      //x+y
      n2[0] = nx+1; n2[1] = ny+1; //torus(nx2, ny2); 
      i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

      //-x-y
      n2[0] = nx-1; n2[1] = ny-1; //torus(nx2, ny2); 
      i2 = j(n2, 0); connectivity[i1][n][2] = i2;

    } else {

      //-x+y
      n2[0] = nx-1; n2[1] = ny+1; //torus(nx2, ny2); 
      i2 = j(n2, 0); links[i1][n][1] = connectivity[i1][n][1] = i2;

      //x-y
      n2[0] = nx+1; n2[1] = ny-1; //torus(nx2, ny2); 
      i2 = j(n2, 0); connectivity[i1][n][2] = i2;
    }
  }
}

/*
void lattice::adjacency(int n)
{
  int **m = createint(LD);

  for(int i=0; i<LD; i++)
    for(int j=0; j<LD; j++)
      m[i][j] = 0;

  for(int i=0; i<LD; i++)
  {
    for(int j=1; j<=links[i][n][0]; j++)
      m[i][links[i][n][j]] = 1;
  }

  cout << "Adjacency matrix of level " << n << ":\n";
  write_m(m, LD);

  destroy(m, LD);
}
*/

void lattice::print() //write out link and connectivity matrices
{
  int Qtmp = SUBL;

  cout << "links:\n";
  for(int i=0; i<Qtmp*LD; i++)
  {
    cout << "i1: " << i;
    for(int n=0; n<nbr; n++) {
      cout << " (" << n << ")";
      for(int j=1; j<=links[i][n][0]; j++) cout << " " << links[i][n][j];
      cout << "; ";
    }
    cout << "\n";
  }
  cout << "\nconnectivity:\n";
  for(int i=0; i<Qtmp*LD; i++)
  {
    cout << "i1: " << i;
    for(int n=0; n<nbr; n++) {
      cout << " (" << n << ")";
      for(int j=1; j<=connectivity[i][n][0]; j++) cout << " " << connectivity[i][n][j];
      cout << "; ";
    }
    cout << "\n";
  }
}

//set a simple periodic chain
void lattice::set_chain()
{
  if( DIM != 1 )
  {
    cout << "lattice::set_chain(): need to compile with DIM=1." << endl;
    exit(-1);
  }
  if( SUBL != 1 )
  {
    cout << "lattice::set_chain(): cannot set chain; number of sites in unit cell Q is different from 1 "<< endl;
    exit(-1);
  }

  desc = "chain";
  cout << "Setting lattice to " << desc << ".\n";

  //first neighbor
  int n=0;
  for(int j1=0; j1<L; j1++) {
    links[j1][n][0] = 1;
    connectivity[j1][n][0] = 2;
    if( j1<L-1 )
      links[j1][n][1] = connectivity[j1][n][1] = j1+1;
    else
      links[j1][n][1] = connectivity[j1][n][1] = 0;
    if( j1>0 )
      connectivity[j1][n][2] = j1-1;
    else
      connectivity[j1][n][2] = L-1;
  }

  //second neighbor
  n=1;
  for(int j1=0; j1<L; j1++) {
    links[j1][n][0] = 1;
    connectivity[j1][n][0] = 2;
    if( j1<L-2 ) {
      links[j1][n][1] = connectivity[j1][n][1] = j1+2;
    } else {
      if( j1==L-2 )
        links[j1][n][1] = connectivity[j1][n][1] = 0;
      else
        links[j1][n][1] = connectivity[j1][n][1] = 1;
    }
    if( j1>1 )
      connectivity[j1][n][2] = j1-2;
    else
      if( j1==0 )
        connectivity[j1][n][2] = L-2;
      else
        connectivity[j1][n][2] = L-1;
  }

  //third neighbor
  n=2;
  for(int j1=0; j1<L; j1++) {
    links[j1][n][0] = 1;
    connectivity[j1][n][0] = 2;
    if( j1<L-3 ) {
      links[j1][n][1] = connectivity[j1][n][1] = j1+3;
    } else {
      if( j1==L-3 )
        links[j1][n][1] = connectivity[j1][n][1] = 0;
      else if( j1==L-2 )
        links[j1][n][1] = connectivity[j1][n][1] = 1;
      else
        links[j1][n][1] = connectivity[j1][n][1] = 2;
    }
    if( j1>2 ) {
      connectivity[j1][n][2] = j1-3;
    } else {
      if( j1==0 )
        connectivity[j1][n][2] = L-3;
      else if( j1==1 )
        connectivity[j1][n][2] = L-2;
      else
        connectivity[j1][n][2] = L-1;
    }
  }
}

