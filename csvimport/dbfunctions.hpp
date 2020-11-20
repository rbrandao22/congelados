#ifndef DBFUNCTIONSHEADERDEF
#define DBFUNCTIONSHEADERDEF

#include <pqxx/pqxx>
#include <string>

void csvimport(pqxx::work& W, std::string datafile, std::string table);

#endif
