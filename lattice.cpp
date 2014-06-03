#include "lattice.h"
#include "helperfunctions.h"

lattice::lattice(int L0, int Q0)
{
  this->L = L0; this->Q = Q0;
  L2 = L*L; N = L2*Q;

  cout << "Creating lattice with L=" << L << " and " << "Q=" << Q << "." << endl;

  //number of neighbors we want to store in the links/connectivity matrices
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
}

// mapping from cartesian coordinates to linear coordinates
int lattice::j(int n0, int m0, int q0)
{
#ifdef SUBL
    return Q*(n0 + m0*L) + q0;
#else
    return n0 + m0*L;
#endif
}

std::string lattice::get_desc() {return desc;}

// inverse mapping from linear to cartesian coordinates
void lattice::getnq(int &nx, int &ny, int &q, int jj)
{
  //cout << "L=" << L << endl;
  //cout << "this->L=" << this->L << endl;

  if(jj<0) {
    cout << "ERROR: undefined mod...\n";
    exit(-1);
  }
#if SUBL
  q = jj % Q;
  int j0 = jj / Q;
  nx = j0 % L;
  ny = j0 / L;
#else
  if( jj==0 )
  {
    nx = 0; ny = 0;
    return;
  }
  nx = jj % L;
  ny = jj / L;
#endif
}

//restrain the values to the LxL using periodc boundary conditions
void lattice::torus(int &nx, int &ny)
{
  if(nx<0 ) nx += L;
  if(nx>=L) nx -= L;

  if(ny<0 ) ny += L;
  if(ny>=L) ny -= L;
}

//sets the connectivity and links matrices to a kagome (currently only nearest neighbor!)
void lattice::set_kagome()
{
  int n;
  int nx2, ny2, q2, i2;
  int i1;

  //checking for consistency
  desc = "kagome";
  cout << "Setting lattice to " << desc << ".\n";

  if( Q!=3 )
  {
    cout << "lattice::set_kagome():: cannot set kagome; number of sites in the unit cell Q is different from 3 "<<endl;
    exit(-1);
  }

  //loop over all UNIT CELLS of the lattice (cell: three sites of the hexagon)
  
  for(int nx=0; nx<L; nx++) // x-coordinate
  {
    for(int ny=0; ny<L; ny++) // y-coordinate
    {
      for(int q=0; q<Q; q++) // loop over positions in this unit cell
      {
        // get the linear index of the current site 
        i1 = j(nx, ny, q);

        //setting nearest neighbors: links =  distinct links in the direction 0, 60 and 120 degree; connectivity = all links
        n = 0;
        links[i1][n][0] = 2; //number of nn for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        { //ordering: straight lignes, then one rotation to the right
          // +y
          nx2 = nx; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // -y
          nx2 = nx; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +x+y
          nx2 = nx+1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          
          // -x-y
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
/*
          // +x
          nx2 = nx; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +y
          nx2 = nx; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x
          nx2 = nx-1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // -y
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;  
*/
        }  

        if( q==1 ) // second site in the unit cell
        {
          // -x
          nx2 = nx; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x
          nx2 = nx+1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +y
          nx2 = nx; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;

          // -y
          nx2 = nx; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;

/*
          // +x 
          nx2 = nx+1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -x+y
          nx2 = nx; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x
          nx2 = nx; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // +x-y
          nx2 = nx+1; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;  
*/
        }

        if( q==2 ) // third site in the unit cell
        {
          // -x-y
          nx2 = nx-1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x+y
          nx2 = nx; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -x
          nx2 = nx-1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;

          // +x
          nx2 = nx; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;

/*
          // +y
          nx2 = nx; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -x+y
          nx2 = nx-1; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -y
          nx2 = nx; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // +x-y
          nx2 = nx; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;  
*/
        }

        // second neighbors
        n = 1;
        links[i1][n][0] = 2;//number of 2nd for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        { //ordering: straight lines, then one rotation to the left

          // -x+y
          nx2 = nx; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x-y
          nx2 = nx+1; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -2x-y
          nx2 = nx-1; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;

          // +2x+y
          nx2 = nx+1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;

/*
          // 2x - y
          nx2 = nx+1; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -x + 2y
          nx2 = nx-1; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -2x + y
          nx2 = nx-1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // x - 2y
          nx2 = nx; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
*/
        }

        if( q==1 ) // second site in the unit cell
        {
          // -2x-y
          nx2 = nx-1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x+y
          nx2 = nx+1; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // -x-2y
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;

          // +x+2y
          nx2 = nx+1; ny2 = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;

/*
          // x + y
          nx2 = nx+1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -x + 2y
          nx2 = nx; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x - y
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // x - 2y
          nx2 = nx+1; ny2 = ny-1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
*/
        }

        if( q==2 ) // third site in the unit cell
        {
          // -x-2y
          nx2 = nx-1; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +x+2y
          nx2 = nx; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;

          // +x-y
          nx2 = nx; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;

          // -x+y
          nx2 = nx-1; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
/*
          // 2x - y
          nx2 = nx+1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // x + y
          nx2 = nx+1; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -2x + y
          nx2 = nx-1; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // -x - y
          nx2 = nx-1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
*/
        }

/* deactivatind non-diagonal third neighbor for now, because we are not using it for now
        // third neighbors (non-diagonal)
        n = 2;
        links[i1][n][0] = 2;  //number of 2nd for a given site
        connectivity[i1][n][0] = 4;

        if( q==0 ) // first site in the unit cell
        {
          // +x,0
          nx2 = nx+1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +x+y,0
          nx2 = nx+1; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x,0
          nx2 = nx-1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // -x-y,0
          nx2 = nx-1; ny2 = ny-1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==1 ) // second site in the unit cell
        {
          // +x,1
          nx2 = nx+1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +y,1
          nx2 = nx; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x,1
          nx2 = nx-1; ny2 = ny; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // -y,1
          nx2 = nx; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
        }

        if( q==2 ) // third site in the unit cell
        {
          // +x+y,2
          nx2 = nx+1; ny2 = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // +y,2
          nx2 = nx; ny2 = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][2] = connectivity[i1][n][2] = i2;
          // -x-y,2
          nx2 = nx-1; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][3] = i2;
          // -y,2
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][4] = i2;
        }
*/
        // diagonal neighbors (inside the hexagon)
        n = 2;
        links[i1][n][0] = 1;//number of 2nd for a given site
        connectivity[i1][n][0] = 2;

        if( q==0 ) // first site in the unit cell
        {
          // -2x
          nx2 = nx-1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x
          nx2 = nx+1; ny2 = ny; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;
/*
          // 2x - 2y
          nx2 = nx+1; ny2 = ny-1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -2x + 2y
          nx2 = nx-1; ny2 = ny+1; q2 = 0; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;
*/
        }

        if( q==1 ) // second site in the unit cell
        {
          // -2x-2y
          nx2 = nx-1; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2x+2y
          nx2 = nx+1; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;
/*
          // 2y
	  x2 = nx; ny2 = ny+1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -2y
          nx2 = nx; ny2 = ny-1; q2 = 1; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;
*/
        }

        if( q==2 ) // third site in the unit cell
        {
          // -2y
          nx2 = nx; ny2 = ny-1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;

          // +2y
          nx2 = nx; ny2 = ny+1; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;

/*
          // 2x
          nx2 = nx+1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); links[i1][n][1] = connectivity[i1][n][1] = i2;
          // -2x
          nx2 = nx-1; ny2 = ny; q2 = 2; torus(nx2, ny2);
          i2 = j(nx2, ny2, q2); connectivity[i1][n][2] = i2;
*/
        }

      } //loop over sites in unit cell
    } //loop over sites (x)
  } //loop over sites (y)
}

//sets the connectivity and links matrices to a triangular lattice (currently up to third neighbor, and four-site plaquette term)
void lattice::set_triangular()
{
  int n;
  int nx2, ny2, i2;
  int nx, ny, q;

  desc = "triangular";
  cout << "Setting lattice to " << desc << ".\n";

  for(int i1=0; i1<L2; i1++)
  {
    getnq(nx, ny, q, i1);

    //nearest neighbor
    n=0;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //x
    nx2 = nx+1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //x+y
    nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //y
    nx2 = nx; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-x
    nx2 = nx-1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //-x-y
    //y direction
    nx2 = nx-1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][5] = i2;

    //-y
    nx2 = nx; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][6] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //2x+y
    nx2 = nx+2; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //x+2y
    nx2 = nx+1; ny2 = ny+2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x+y
    nx2 = nx-1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-2x-y
    nx2 = nx-2; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //-x-2y
    nx2 = nx-1; ny2 = ny-2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][5] = i2;

    //x-y
    nx2 = nx+1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][6] = i2;

    //3rd neighbor
    n=2;
    links[i1][n][0] = 3;
    connectivity[i1][n][0] = 6;

    //2x
    nx2 = nx+2; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //2x+2y
    nx2 = nx+2; ny2 = ny+2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //2y
    nx2 = nx; ny2 = ny+2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][3] = connectivity[i1][n][3] = i2;

    //-2x
    nx2 = nx-2; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //-2x-2y
    nx2 = nx-2; ny2 = ny-2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][5] = i2;

    //-2y
    nx2 = nx; ny2 = ny-2; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][6] = i2;

    //4-site plaquette term
    n=0;
    plaquette4[i1][n][0] = 3;
    nx2 = nx+1; ny2 = ny  ; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][1] = i2;
    nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][2] = i2;
    nx2 = nx  ; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][3] = i2;

    nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][4] = i2;
    nx2 = nx+1; ny2 = ny+2; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][5] = i2;
    nx2 = nx  ; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][6] = i2;

    nx2 = nx+1; ny2 = ny  ; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][7] = i2;
    nx2 = nx+2; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][8] = i2;
    nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][9] = i2;
  }
}

//square lattice on a torus
void lattice::set_square()
{
  int n, nx2, ny2, i2;
  int nx, ny, q;

  desc = "square";
  cout << "Setting lattice to " << desc << ".\n";

  for(int i1=0; i1<L2; i1++)
  {
    getnq(nx, ny, q, i1);

    //nearest neighbor
    n=0;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x
    nx2 = nx+1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //y
    nx2 = nx; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x
    nx2 = nx-1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][3] = i2;

    //-y
    nx2 = nx; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x+y
    nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //-x+y
    nx2 = nx-1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x-y
    nx2 = nx-1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][3] = i2;

    //x-y
    nx2 = nx+1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //4-site plaquette term
    n=0;
    plaquette4[i1][n][0] = 1;
    nx2 = nx+1; ny2 = ny2  ; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][1] = i2;
    nx2 = nx+1; ny2 = ny2+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][2] = i2;
    nx2 = nx  ; ny2 = ny2+1; torus(nx2, ny2); i2 = j(nx2,ny2,0); plaquette4[i1][n][3] = i2;
  }
}

void lattice::set_checkerboard()
{
  int n, nx2, ny2, i2;

  if(L%2!=0) {
    cout << "Error: cannot set checkerboard... (L=" << L << ")\n";
    exit(-1);
  }

  desc = "checkerboard";
  cout << "Setting lattice to " << desc << ".\n";

  int nx, ny, q;

  for(int i1=0; i1<L2; i1++)
  {
    getnq(nx, ny, q, i1);

    //nearest neighbor
    n=0;
    links[i1][n][0] = 2;
    connectivity[i1][n][0] = 4;

    //x
    nx2 = nx+1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][1] = connectivity[i1][n][1] = i2;

    //y
    nx2 = nx; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    links[i1][n][2] = connectivity[i1][n][2] = i2;

    //-x
    nx2 = nx-1; ny2 = ny; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][3] = i2;

    //-y
    nx2 = nx; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
    connectivity[i1][n][4] = i2;

    //2nd neighbor
    n=1;
    links[i1][n][0] = 1;
    connectivity[i1][n][0] = 2;

    if( ((nx+ny)%2)==0 )
    { //right

      //x+y
      nx2 = nx+1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
      links[i1][n][1] = connectivity[i1][n][1] = i2;

      //-x-y
      nx2 = nx-1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
      connectivity[i1][n][2] = i2;

    } else {

      //-x+y
      nx2 = nx-1; ny2 = ny+1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
      links[i1][n][1] = connectivity[i1][n][1] = i2;

      //x-y
      nx2 = nx+1; ny2 = ny-1; torus(nx2, ny2); i2 = j(nx2,ny2,0);
      connectivity[i1][n][2] = i2;
    }
  }
}

/*
void lattice::adjacency(int n)
{
  int **m = createint(L2);

  for(int i=0; i<L2; i++)
    for(int j=0; j<L2; j++)
      m[i][j] = 0;

  for(int i=0; i<L2; i++)
  {
    for(int j=1; j<=links[i][n][0]; j++)
      m[i][links[i][n][j]] = 1;
  }

  cout << "Adjacency matrix of level " << n << ":\n";
  write_m(m, L2);

  destroy(m, L2);
}
*/

void lattice::print() //write out link and connectivity matrices
{

#ifdef SUBL
  int Qtmp = Q;
#else
  int Qtmp = 1;
#endif

  cout << "links:\n";
  for(int i=0; i<Qtmp*L2; i++)
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
  for(int i=0; i<Qtmp*L2; i++)
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

