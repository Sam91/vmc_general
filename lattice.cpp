#include "lattice.h"
#include "helperfunctions.h"

lattice::lattice(int L0, int Q0)
{
  this->L = L0; this->Q = Q0; this->L2 = L0*L0;
  //nx = ny = q = 0;
  neighbors = new int[L2];

  nbr=3;

  /* The site relations are defined in the following way
   * plaquette[site index j1][type of plaquette][number of plaquettes on this site, j2, j3, j4, ...]
  */

  connectivity = new (int(**[this->L2])); //lets take first and second neighbor
  links = new (int(**[this->L2]));
  plaquette4 = new (int(**[this->L2])); //four-site plaquettes

  for(int i=0; i<L2; i++)
  {
    connectivity[i] = createint(nbr,10);
    links[i] = createint(nbr,10);
    plaquette4[i] = createint(nbr,20);
  }

}

lattice::~lattice()
{
  delete[] neighbors;

  for(int i=0; i<L2; i++)
  {
    destroy(connectivity[i], nbr);
    destroy(links[i], nbr);
    destroy(plaquette4[i], nbr);
  }
  delete[] connectivity;
  delete[] links;
  delete[] plaquette4;
}

int lattice::j(int n0, int m0, int q0)
{
#ifdef SUBL
    return n0 + m0*L + q0*L2;
#else
    return n0 + m0*L;
#endif
}

std::string lattice::get_desc() {return desc;}

void lattice::getnq(int &nx, int &ny, int &q, int jj)
{
  //cout << "L=" << L << endl;
  //cout << "this->L=" << this->L << endl;

  if(jj<0) {
    cout << "ERROR: undefined mod...\n";
    exit(-1);
  }
#if SUBL
  q = jj % L2;
  int j0 = jj - L2*q;
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
void lattice::torus( int &nx, int &ny)
{
  if(nx<0 ) nx += L;
  if(nx>=L) nx -= L;

  if(ny<0 ) ny += L;
  if(ny>=L) ny -= L;
}
/*
void lattice::torus()
{
  if(nx<0) nx += L;
  nx = (nx%L);
  if(ny<0) ny += L;
  ny = (ny%L);
}*/

//sets the connectivity and links matrices to a triangular one (currently only next neighbor!)
//triangular lattice on a torus
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

//finds all n-th neighbors to site i1 and writes them to table "neighbors".
//Here, the connectivity of the lattice is defined
void lattice::findn(int i1, int n)
{
  switch( n )
  {
    case 1:
      neighbors[0] = 4;

      //x-direction
      if( i1%L == L-1 ) neighbors[1] = i1-L+1;
      else neighbors[1] = i1+1;

      if( i1%L == 0 ) neighbors[2] = i1+L-1;
      else neighbors[2] = i1-1;

      //y-direction
      if( i1>=L2-L ) neighbors[3] = i1-L2+L;
      else neighbors[3] = i1+L;

      if( i1<L ) neighbors[4] = i1+L2-L;
      else neighbors[4] = i1-L;

      break;
    case 2:
      neighbors[0] = 4;

      //xy-next nearest
      if( i1==L2-1 ) neighbors[1] = 0;
      else if( i1>=L2-L ) neighbors[1] = i1-L2+L+1;
      else if( i1%L==L-1 ) neighbors[1] = i1+1;
      else neighbors[1] = i1+L+1;

      if( i1==0 ) neighbors[2] = L2-1;
      else if( i1<=L ) neighbors[2] = i1+L2-L-1;
      else if( i1%L==0 ) neighbors[2] = i1-1;
      else neighbors[2] = i1-L-1;

      //-xy-next nearest
      if( i1==L2-L ) neighbors[3] = L-1;
      else if( i1>L2-L ) neighbors[3] = i1-L2+L-1;
      else if( i1%L == 0 ) neighbors[3] = i1+2*L-1;
      else neighbors[3] = i1+L-1;

      if( i1==L-1 ) neighbors[4] = L2-L;
      else if( i1<L ) neighbors[4] = i1+L2-L+1;
      else if( i1%L == L-1 ) neighbors[4] = i1-2*L+1;
      else neighbors[4] = i1-L+1;

      break;
    default:
      neighbors[0] = 0;
      //std::cout << "Error\n";
  }
}

void lattice::findn0(int i1, int n)
{
  switch( n )
  {
    case 1:
      neighbors[0] = 2;

      //x-direction
      if( i1%L == L-1 ) neighbors[1] = i1-L+1;
      else neighbors[1] = i1+1;

      //y-direction
      if( i1>=L2-L ) neighbors[2] = i1-L2+L;
      else neighbors[2] = i1+L;

      break;

    case 2:
      neighbors[0] = 2;

      //xy-next nearest
      if( i1==L2-1 ) neighbors[1] = 0;
      else if( i1>=L2-L ) neighbors[1] = i1-L2+L+1;
      else if( i1%L==L-1 ) neighbors[1] = i1+1;
      else neighbors[1] = i1+L+1;

      //-xy-next nearest
      if( i1==L2-L ) neighbors[2] = L-1;
      else if( i1>L2-L ) neighbors[2] = i1-L2+L-1;
      else if( i1%L == 0 ) neighbors[2] = i1+2*L-1;
      else neighbors[2] = i1+L-1;

      break;

    default:
      neighbors[0] = 0;
      //std::cout << "Error\n";
  }
}

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

void lattice::print() //write out link and connectivity matrices
{
  cout << "links:\n";
  for(int i=0; i<L2; i++)
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
  for(int i=0; i<L2; i++)
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
