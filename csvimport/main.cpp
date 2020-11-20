#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>

#include "dbfunctions.hpp"

using std::string;
using std::vector;


int main(int argc, char* argv[])
{
    try {
        pqxx::connection C("dbname = congelados user = postgres password = "\
			   "passwd hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() <<\
	      std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }

        // csv import argv[1]
        if (argc > 1 && (std::strcmp(argv[1], "dicionario") == 0 ||\
			 std::strcmp(argv[1], "populacao") == 0 ||\
			 std::strcmp(argv[1], "renda") == 0 ||\
			 std::strcmp(argv[1], "idade") == 0 ||\
			 std::strcmp(argv[1], "lanches_shares") == 0 ||\
			 std::strcmp(argv[1], "lanches_precos") == 0)) {
            pqxx::work W(C);
	    string table = argv[1];
            string datafile = "/var/lib/postgresql/data/pgdata/databases/" \
	      "congelados/" + table + ".csv";
            csvimport(W, datafile, table);
            W.commit();
	} else {
	  std::cout << "Please provide valid table name as argument" <<\
	    std::endl;
	}

        std::cout << "database operation successfull" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
