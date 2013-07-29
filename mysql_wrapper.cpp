#include "mysql_wrapper.h"
#include "vmc.h"
#include <sstream>
#include <cstdlib>
#include <time.h>

//The DB connection parameters are specified here
mysql_wrapper::mysql_wrapper()
{

/*  serverString = "tcp://10.1.1.1:3306";
  user = "vmc";
  passwd = "vmc";
  db = "vmc";
*/

  #include "mysql_user.txt"

  driver = getdriver();
  con = getconnection();
  //stmt=0; con=0; res=0;
}

mysql_wrapper::~mysql_wrapper()
{
  if(con!=0) {delete con; con=0;}
}

//tries to get a driver several times
sql::Driver* mysql_wrapper::getdriver() {
  sql::Driver* driver;

  for(int p=0; p<10; p++) {
    try {
      driver = get_driver_instance();
    } catch (sql::SQLException &e) {
      std::cout << "# ERR: SQLException in " << __FILE__;
      std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
      std::cout << "# ERR: " << e.what();
      std::cout << " (MySQL error code: " << e.getErrorCode();
      std::cout << ", SQLState: " << e.getSQLState() << " )\n";

      sleep(1);
      continue;
    }
    std::cout << __FUNCTION__ << ": we've got a mysql driver\n";
    return driver;
  }
  std::cout << __FUNCTION__ << ": no mysql driver gotten\n";
  exit(-1);
}

sql::Connection* mysql_wrapper::getconnection() {
  sql::Connection* connection;

  for(int p=0; p<10; p++) {
    try {
      connection = driver->connect(serverString, user, passwd);
      connection->setSchema( db );
    } catch (sql::SQLException &e) {
      std::cout << "# ERR: SQLException in " << __FILE__;
      std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
      std::cout << "# ERR: " << e.what();
      std::cout << " (MySQL error code: " << e.getErrorCode();
      std::cout << ", SQLState: " << e.getSQLState() << " )\n";

      sleep(1);
      continue;
    }
    std::cout << __FUNCTION__ << ": we've got a mysql connection\n";
    return connection;
  }
  std::cout << __FUNCTION__ << ": no mysql connection gotten\n";
  exit(-1);
}


//insert to s-wave table (we probably need a table for each state)
int mysql_wrapper::insert_s(int* sys, double* values)
{
  std::ostringstream os;

  //build up the query
  os << "INSERT INTO singletsc VALUES (''";
  for(int i=0; i<2; i++)
    os << "," << sys[i];
  for(int i=0; i<19; i++)
    os << "," << values[i];
  os << ")";
  std::cout << os.str() << "\n";

  //driver = getdriver();
  //con = getconnection();

  try {

    stmt = con->createStatement();
    stmt->executeUpdate( os.str() );

//    if(stmt!=0) {delete stmt; stmt=0;}
//    if(con!=0) {delete con; con=0;}
    delete stmt;

  } catch (sql::SQLException &e) {
    std::cout << "# ERR: SQLException in " << __FILE__;
    std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
    std::cout << "# ERR: " << e.what();
    std::cout << " (MySQL error code: " << e.getErrorCode();
    std::cout << ", SQLState: " << e.getSQLState() << " )\n";
    delete stmt;
    
    return -1;
  }

  return 0;
}

//retrieve from s-wave table for given sys and variational parameters values[0], values[1]
int mysql_wrapper::get_s(int* sys, double* values)
{
  int ret;
  std::ostringstream os;

  //build up the query
  os << "SELECT * FROM singletsc WHERE sites=" << sys[0] << " AND Nz=" << sys[1] << " AND mux=" << values[0] << " AND dds=" << values[1] << " AND ddd=" << values[2];
  std::cout << os.str() << "\n";

  //driver = getdriver();
  //con = getconnection();

  try {
//    driver = get_driver_instance();
//    con = driver->connect(serverString, user, passwd);

    stmt = con->createStatement();
    res = stmt->executeQuery( os.str() );

    //number of real variational parameters
    int nq = 3;

    if( res->next() ) {
      values[nq   ] = res->getDouble("SS");
      values[nq+ 1] = res->getDouble("dSS");
      values[nq+ 2] = res->getDouble("SS_mux");
      values[nq+ 3] = res->getDouble("dSS_mux");
      values[nq+ 4] = res->getDouble("SS_dds");
      values[nq+ 5] = res->getDouble("dSS_dds");
      values[nq+ 6] = res->getDouble("SS2");
      values[nq+ 7] = res->getDouble("dSS2");
      values[nq+ 8] = res->getDouble("SS2_mux");
      values[nq+ 9] = res->getDouble("dSS2_mux");
      values[nq+10] = res->getDouble("SS2_dds");
      values[nq+11] = res->getDouble("dSS2_dds");
      values[nq+12] = res->getDouble("SS_ddd");
      values[nq+13] = res->getDouble("dSS_ddd");
      values[nq+14] = res->getDouble("SS2_ddd");
      values[nq+15] = res->getDouble("dSS2_ddd");
      ret = 1;
    } else ret = 0;

    delete res; delete stmt;

  } catch (sql::SQLException &e) {
    std::cout << "# ERR: SQLException in " << __FILE__;
    std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
    std::cout << "# ERR: " << e.what();
    std::cout << " (MySQL error code: " << e.getErrorCode();
    std::cout << ", SQLState: " << e.getSQLState() << " )\n";

    ret = -1;
  }
//  if(stmt!=0) {delete stmt; stmt=0;}
//  if(con!=0) {delete con; con=0;}
//  if(res!=0) {delete res; res=0;}

  return ret;
}

//save and retrieve simple values (without gradient terms)
//insert to evenz1 table (unpaired z-fermions and even pairing: s, ext-s, d)
//parameters: sites, Nz, mux, dds, ddse, ddd
/*
int mysql_wrapper::insert_evenz1(parameters &parms, double* values)
{
  std::ostringstream os;

  os << "INSERT INTO evenz1 (id, sites, Nz, ph, mux, dds, ddse, ddd, ddfxy, ddpxy, muz, ddfz, ddpz, p, k, ord3, ord3a, ord3b, ";
  os << "AUx, AUy, AUz, AVx, AVy, AVz, BUx, BUy, BUz, BVx, BVy, BVz, CUx, CUy, CUz, CVx, CVy, CVz, comment, ";
  os << "SS, dSS, SS2, dSS2, R3, dR3, SSn, dSSn, SS2n, dSS2n, ";
  os << "ANx, dANx, ANz, dANz, ANy, dANy, BNx, dBNx, BNz, dBNz, BNy, dBNy, CNx, dCNx, CNz, dCNz, CNy, dCNy, ";
  os << "Qxy, dQxy, Qxz, dQxz, Q, dQ, gxy, dgxy, gxz, dgxz, ";
  os << "Qxyn, dQxyn, Qxzn, dQxzn, Qn, dQn, gxyn, dgxyn, gxzn, dgxzn ";
  os << ") VALUES (''";

  //system parameters
  os << "," << parms.L2 << "," << parms.Nz << "," << parms.ph;

  //variational parameters (do some rounding)
  os << ",round(" << parms.mux << ",2)" << ",round(" << parms.dds  << ",2)" << ",round(" << parms.ddse << ",2)" << ",round(" << parms.ddd << ",2)";
  os << ",round(" << parms.ddfxy << ",2)" << ",round(" << parms.ddpxy << ",2)";
  os << ",round(" << parms.muz << ",2)" << ",round(" << parms.ddfz << ",2)" << ",round(" << parms.ddpz << ",2)";
  os << ",round(" << parms.p << ",2)";
  os << ",round(" << parms.k << ",2)";
  os << ",round(" << parms.ord3 << ",2)";
  os << ",round(" << parms.ord3a << ",2)";
  os << ",round(" << parms.ord3b << ",2)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.A.getU(i) << ",6)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.A.getV(i) << ",6)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.B.getU(i) << ",6)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.B.getV(i) << ",6)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.C.getU(i) << ",6)";
  for(int i=0; i<3; i++) os << ",round(" << parms.mf.C.getV(i) << ",6)";
  if( parms.comment.size()>0 ) os << ", '" << parms.comment << "'";
  else os << ",''";

  //values
  for(int i=0; i<10; i++)
    os << "," << values[i];

  for(int i=2*91; i<2*100; i++)
    os << "," << values[i];

//  for(int i=2*65; i<2*70; i++)
  for(int i=2*136; i<2*146; i++)
    os << "," << values[i];

  os << ")";
  std::cout << os.str() << "\n";

  try
  {
    stmt = con->createStatement();
    stmt->executeUpdate( os.str() );

    delete stmt;

  } catch (sql::SQLException &e) {
    std::cout << "# ERR: SQLException in " << __FILE__;
    std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
    std::cout << "# ERR: " << e.what();
    std::cout << " (MySQL error code: " << e.getErrorCode();
    std::cout << ", SQLState: " << e.getSQLState() << " )\n";

    return -1;
  }

  return 0;
}

//int mysql_wrapper::get_evenz1(int* sys, double* values, bool withK3)
int mysql_wrapper::get_evenz1(parameters &parms, double* values, bool withQ)
{
  int ret;
  std::ostringstream os;

  //build up the query
  os << "SELECT * FROM evenz1 WHERE sites=" << parms.L2 << " AND Nz=" << parms.Nz << " AND ph=" << parms.ph;
  os << " AND mux=round(" << parms.mux << ",2) AND dds=round(" << parms.dds  << ",2) AND ddse=round(" << parms.ddse << ",2) AND ddd=round(" << parms.ddd << ",2)";
  os << " AND ddfxy=round(" << parms.ddfxy << ",2) AND ddpxy=round(" << parms.ddpxy << ",2)";
  os << " AND muz=round(" << parms.muz << ",2) AND ddfz=round(" << parms.ddfz << ",2) AND ddpz=round(" << parms.ddpz << ",2)";
  os << " AND round(AUx,5)=round(" << parms.mf.A.getU(0) << ",5) AND round(AUy,5)=round(" << parms.mf.A.getU(1) << ",5) AND round(AUz,5)=round(" << parms.mf.A.getU(2) << ",5)";
  os << " AND round(AVx,5)=round(" << parms.mf.A.getV(0) << ",5) AND round(AVy,5)=round(" << parms.mf.A.getV(1) << ",5) AND round(AVz,5)=round(" << parms.mf.A.getV(2) << ",5)";
  os << " AND round(BUx,5)=round(" << parms.mf.B.getU(0) << ",5) AND round(BUy,5)=round(" << parms.mf.B.getU(1) << ",5) AND round(BUz,5)=round(" << parms.mf.B.getU(2) << ",5)";
  os << " AND round(BVx,5)=round(" << parms.mf.B.getV(0) << ",5) AND round(BVy,5)=round(" << parms.mf.B.getV(1) << ",5) AND round(BVz,5)=round(" << parms.mf.B.getV(2) << ",5)";
  os << " AND round(CUx,5)=round(" << parms.mf.C.getU(0) << ",5) AND round(CUy,5)=round(" << parms.mf.C.getU(1) << ",5) AND round(CUz,5)=round(" << parms.mf.C.getU(2) << ",5)";
  os << " AND round(CVx,5)=round(" << parms.mf.C.getV(0) << ",5) AND round(CVy,5)=round(" << parms.mf.C.getV(1) << ",5) AND round(CVz,5)=round(" << parms.mf.C.getV(2) << ",5)";
  if( withQ ) os << " AND Q IS NOT NULL";

  std::cout << os.str() << "\n";

  try {
    stmt = con->createStatement();
    res = stmt->executeQuery( os.str() );

    //number of real variational parameters
    int nq = 0;

    if( res->next() ) {
      values[nq   ] = res->getDouble("SS"  );
      values[nq+ 1] = res->getDouble("dSS" );
      values[nq+ 2] = res->getDouble("SS2" );
      values[nq+ 3] = res->getDouble("dSS2");
      values[nq+ 4] = res->getDouble("R3");
      values[nq+ 5] = res->getDouble("dR3");
      if( withQ ) {
        //values[nq+ 4] = res->getDouble("SS3" );
        //values[nq+ 5] = res->getDouble("dSS3");
        values[nq+ 6] = res->getDouble("Q");
        values[nq+ 7] = res->getDouble("dQ");
      }
      ret = 1;
    } else ret = 0;

    delete res; delete stmt;

  } catch (sql::SQLException &e) {
    std::cout << "# ERR: SQLException in " << __FILE__;
    std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
    std::cout << "# ERR: " << e.what();
    std::cout << " (MySQL error code: " << e.getErrorCode();
    std::cout << ", SQLState: " << e.getSQLState() << " )\n";

    ret = -1;
  }
  return ret;
}
*/

int mysql_wrapper::insert_qry(std::string s)
{
  std::cout << "Wrapper: '" << s << "'\n";

  try
  {
    stmt = con->createStatement();
    stmt->executeUpdate( s );

    delete stmt;

  } catch (sql::SQLException &e)
  {
    std::cout << "# ERR: SQLException in " << __FILE__;
    std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << "\n";
    std::cout << "# ERR: " << e.what();
    std::cout << " (MySQL error code: " << e.getErrorCode();
    std::cout << ", SQLState: " << e.getSQLState() << " )\n";
    delete stmt;

    return -1;
  }

  return 0;
}

