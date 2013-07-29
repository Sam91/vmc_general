#ifndef MYSQL_WRAPPER_H
#define MYSQL_WRAPPER_H

//#include <mysql_connection.h>
#include <string.h>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>

#include "vmc.h"

class mysql_wrapper{

public:
  mysql_wrapper();
  ~mysql_wrapper();

  int insert_s(int*, double *);
  int get_s(int*, double *);

  //int insert_evenz1(parameters &p, double *values);
  //int get_evenz1(parameters &p, double *values, bool);

  int insert_qry(std::string);

  //gets a db driver with some time delay
  sql::Driver* getdriver();
  sql::Connection* getconnection();

private:
  std::string serverString;// = "tcp://172.16.200.210:3306";
  std::string db;// = "vmc_triangular";
  std::string user;// = "vmc";
  std::string passwd;// = "vmc";

  sql::Driver *driver;
  sql::Connection *con;
  sql::Statement *stmt;
  sql::ResultSet *res; 

};

#endif
